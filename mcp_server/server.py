#!/usr/bin/env python3
"""
MCP Server for nf-treeseq Nextflow Pipeline Management

This server implements the Model Context Protocol (MCP) to provide AI assistants
with tools for managing Nextflow pipelines.

Protocol: Line-delimited JSON over stdin/stdout (JSON-RPC style)

Client messages:
  {"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"test-client","version":"1.0"}}}
  {"jsonrpc":"2.0","id":2,"method":"tools/list"}
  {"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"listWorkflows","arguments":{}}}

Server responses:
  {"jsonrpc":"2.0","id":1,"result":{"protocolVersion":"2024-11-05","capabilities":{"tools":{}},"serverInfo":{"name":"nf-treeseq-mcp","version":"0.1.0"}}}
  {"jsonrpc":"2.0","id":2,"result":{"tools":[...]}}
  {"jsonrpc":"2.0","id":3,"result":{"content":[{"type":"text","text":"..."}]}}

Tools implemented:
- listWorkflows: scans workflows/ and subworkflows/ for .nf files
- launchRun: spawns a Nextflow run and returns a runId
- monitorRun: checks status of an active run
- listRuns: shows all active/recent runs

Safety features:
- Paths confined to repository root
- Only 'nextflow run' command allowed (no arbitrary shell)
- Run isolation in .runs/ directory
- Process tracking and cleanup
"""
from __future__ import annotations
import sys
import json
import uuid
import time
import subprocess
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional
from datetime import datetime

# Configure logging to stderr (stdout is for JSON-RPC)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    stream=sys.stderr
)
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent
WORKFLOWS_DIR = REPO_ROOT / "workflows"
SUBWORKFLOWS_DIR = REPO_ROOT / "subworkflows"
RUNS_DIR = REPO_ROOT / ".runs"
RUNS_DIR.mkdir(exist_ok=True)

ACTIVE_RUNS: Dict[str, Dict[str, Any]] = {}

SERVER_INFO = {
    "name": "nf-treeseq-mcp",
    "version": "0.1.0"
}

PROTOCOL_VERSION = "2024-11-05"


class MCPError(Exception):
    """Base exception for MCP tool errors"""
    pass


def send_response(request_id: Any, result: Any = None, error: Any = None) -> None:
    """Send a JSON-RPC 2.0 response"""
    response = {"jsonrpc": "2.0", "id": request_id}
    if error:
        response["error"] = error
    else:
        response["result"] = result
    sys.stdout.write(json.dumps(response) + "\n")
    sys.stdout.flush()
    logger.debug(f"Sent response: {response}")


def send_notification(method: str, params: Any) -> None:
    """Send a JSON-RPC 2.0 notification (no id)"""
    notification = {"jsonrpc": "2.0", "method": method, "params": params}
    sys.stdout.write(json.dumps(notification) + "\n")
    sys.stdout.flush()


# ----------------------- TOOL IMPLEMENTATIONS ----------------------- #

def tool_listWorkflows(args: Dict[str, Any]) -> Dict[str, Any]:
    """
    List all available Nextflow workflows and subworkflows in the repository.

    Returns:
        {
            "workflows": [
                {"name": "nf-treeseq", "path": "workflows/nf-treeseq.nf", "kind": "workflow"},
                ...
            ],
            "count": 2
        }
    """
    logger.info("Listing workflows...")
    workflows = []

    # Main workflows
    if WORKFLOWS_DIR.exists():
        for item in WORKFLOWS_DIR.glob("*.nf"):
            workflows.append({
                "name": item.stem,
                "path": str(item.relative_to(REPO_ROOT)),
                "kind": "workflow"
            })

    # Subworkflows (recursive)
    if SUBWORKFLOWS_DIR.exists():
        for item in SUBWORKFLOWS_DIR.rglob("*.nf"):
            workflows.append({
                "name": item.stem,
                "path": str(item.relative_to(REPO_ROOT)),
                "kind": "subworkflow"
            })

    logger.info(f"Found {len(workflows)} workflows/subworkflows")
    return {"workflows": workflows, "count": len(workflows)}


def tool_launchRun(args: Dict[str, Any]) -> Dict[str, Any]:
    """
    Launch a Nextflow pipeline run.

    Args:
        profile: Nextflow profile to use (default: "test")
        entry: Path to workflow file (default: "workflows/nf-treeseq.nf")
        params: Dictionary of parameters to pass to the pipeline
        revision: Optional git revision/tag

    Returns:
        {
            "runId": "uuid",
            "status": "STARTED",
            "pid": 12345,
            "entry": "workflows/nf-treeseq.nf",
            "profile": "test",
            "startTime": "2024-11-14T10:30:00"
        }
    """
    profile = args.get("profile", "test")
    revision = args.get("revision")
    entry = args.get("entry", "workflows/nf-treeseq.nf")
    params: Dict[str, Any] = args.get("params", {})

    logger.info(f"Launching run: entry={entry}, profile={profile}")

    # Validate entry path
    entry_path = (REPO_ROOT / entry).resolve()
    if not entry_path.exists() or not entry_path.is_file():
        raise MCPError(f"Entry workflow not found: {entry}")
    if entry_path.suffix != ".nf":
        raise MCPError("Entry must be a .nf file")

    # Security check: ensure path is within repo
    try:
        entry_path.relative_to(REPO_ROOT)
    except ValueError:
        raise MCPError("Entry path must be within repository")

    run_id = str(uuid.uuid4())
    run_dir = RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    # Serialize params to file
    params_file = run_dir / "params.json"
    with params_file.open("w", encoding="utf-8") as fh:
        json.dump(params, fh, indent=2)

    # Build Nextflow command
    cmd = ["nextflow", "run", str(entry_path), "-work-dir", str(run_dir / "work")]
    if profile:
        cmd.extend(["-profile", profile])
    if revision:
        cmd.extend(["-r", revision])

    # Add params as command line arguments
    for k, v in params.items():
        cmd.append(f"--{k}={v}")

    # Prepare log files
    stdout_log = run_dir / "nextflow.stdout.log"
    stderr_log = run_dir / "nextflow.stderr.log"

    # Start process (non-blocking)
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(REPO_ROOT),
            stdout=stdout_log.open("w"),
            stderr=stderr_log.open("w"),
            text=True,
        )
    except FileNotFoundError:
        raise MCPError("'nextflow' command not found in PATH. Please install Nextflow.")

    start_time = datetime.now().isoformat()

    ACTIVE_RUNS[run_id] = {
        "cmd": cmd,
        "pid": proc.pid,
        "proc": proc,
        "start_time": start_time,
        "status": "RUNNING",
        "run_dir": str(run_dir),
        "stdout_log": str(stdout_log),
        "stderr_log": str(stderr_log),
        "workflow": entry,
        "profile": profile,
        "params": params,
    }

    logger.info(f"Started run {run_id} with PID {proc.pid}")

    return {
        "runId": run_id,
        "status": "STARTED",
        "pid": proc.pid,
        "entry": entry,
        "profile": profile,
        "startTime": start_time,
        "runDir": str(run_dir)
    }


def tool_monitorRun(args: Dict[str, Any]) -> Dict[str, Any]:
    """
    Monitor the status of a running pipeline.

    Args:
        runId: The run ID to monitor
        includeLogs: Whether to include recent log lines (default: False)

    Returns:
        {
            "runId": "uuid",
            "status": "RUNNING" | "COMPLETED" | "FAILED",
            "pid": 12345,
            "startTime": "...",
            "endTime": "..." (if finished),
            "exitCode": 0 (if finished),
            "logs": ["line1", "line2"] (if includeLogs=True)
        }
    """
    run_id = args.get("runId")
    include_logs = args.get("includeLogs", False)

    if not run_id:
        raise MCPError("runId is required")

    if run_id not in ACTIVE_RUNS:
        raise MCPError(f"Run {run_id} not found")

    run_info = ACTIVE_RUNS[run_id]
    proc = run_info.get("proc")

    # Check process status
    if proc:
        poll_result = proc.poll()
        if poll_result is not None:
            # Process has finished
            end_time = datetime.now().isoformat()
            run_info["status"] = "COMPLETED" if poll_result == 0 else "FAILED"
            run_info["end_time"] = end_time
            run_info["exit_code"] = poll_result
            run_info.pop("proc", None)  # Remove process object

    result = {
        "runId": run_id,
        "status": run_info["status"],
        "pid": run_info["pid"],
        "startTime": run_info["start_time"],
        "workflow": run_info["workflow"],
        "profile": run_info["profile"],
    }

    if "end_time" in run_info:
        result["endTime"] = run_info["end_time"]
        result["exitCode"] = run_info["exit_code"]

    # Include recent logs if requested
    if include_logs:
        stdout_log = Path(run_info["stdout_log"])
        if stdout_log.exists():
            with stdout_log.open("r") as f:
                lines = f.readlines()
                result["logs"] = [line.rstrip() for line in lines[-20:]]  # Last 20 lines
        else:
            result["logs"] = []

    logger.info(f"Monitored run {run_id}: status={result['status']}")
    return result


def tool_listRuns(args: Dict[str, Any]) -> Dict[str, Any]:
    """
    List all tracked runs.

    Returns:
        {
            "runs": [
                {"runId": "...", "status": "...", "workflow": "...", ...},
                ...
            ],
            "count": 3
        }
    """
    runs = []
    for run_id, run_info in ACTIVE_RUNS.items():
        # Update status by polling if still running
        proc = run_info.get("proc")
        if proc:
            poll_result = proc.poll()
            if poll_result is not None:
                run_info["status"] = "COMPLETED" if poll_result == 0 else "FAILED"
                run_info["end_time"] = datetime.now().isoformat()
                run_info["exit_code"] = poll_result
                run_info.pop("proc", None)

        runs.append({
            "runId": run_id,
            "status": run_info["status"],
            "workflow": run_info["workflow"],
            "profile": run_info["profile"],
            "startTime": run_info["start_time"],
            "pid": run_info["pid"],
        })

    logger.info(f"Listed {len(runs)} runs")
    return {"runs": runs, "count": len(runs)}


# Tool registry with schemas
TOOLS = {
    "listWorkflows": {
        "handler": tool_listWorkflows,
        "schema": {
            "name": "listWorkflows",
            "description": "List all available Nextflow workflows and subworkflows in the repository",
            "inputSchema": {
                "type": "object",
                "properties": {},
                "required": []
            }
        }
    },
    "launchRun": {
        "handler": tool_launchRun,
        "schema": {
            "name": "launchRun",
            "description": "Launch a Nextflow pipeline run with specified parameters",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "profile": {
                        "type": "string",
                        "description": "Nextflow profile to use (e.g., 'test', 'docker')",
                        "default": "test"
                    },
                    "entry": {
                        "type": "string",
                        "description": "Path to the workflow file",
                        "default": "workflows/nf-treeseq.nf"
                    },
                    "params": {
                        "type": "object",
                        "description": "Parameters to pass to the pipeline"
                    },
                    "revision": {
                        "type": "string",
                        "description": "Git revision or tag to run"
                    }
                },
                "required": []
            }
        }
    },
    "monitorRun": {
        "handler": tool_monitorRun,
        "schema": {
            "name": "monitorRun",
            "description": "Monitor the status of a running pipeline",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "runId": {
                        "type": "string",
                        "description": "The run ID to monitor"
                    },
                    "includeLogs": {
                        "type": "boolean",
                        "description": "Include recent log lines in the response",
                        "default": False
                    }
                },
                "required": ["runId"]
            }
        }
    },
    "listRuns": {
        "handler": tool_listRuns,
        "schema": {
            "name": "listRuns",
            "description": "List all tracked pipeline runs",
            "inputSchema": {
                "type": "object",
                "properties": {},
                "required": []
            }
        }
    }
}


# ----------------------- MESSAGE HANDLERS ----------------------- #

def handle_initialize(request_id: Any, params: Dict[str, Any]) -> None:
    """Handle initialize request"""
    client_info = params.get("clientInfo", {})
    logger.info(f"Initializing connection with client: {client_info.get('name', 'unknown')}")

    result = {
        "protocolVersion": PROTOCOL_VERSION,
        "capabilities": {
            "tools": {}
        },
        "serverInfo": SERVER_INFO
    }
    send_response(request_id, result)


def handle_tools_list(request_id: Any) -> None:
    """Handle tools/list request"""
    tools_list = [tool_info["schema"] for tool_info in TOOLS.values()]
    send_response(request_id, {"tools": tools_list})


def handle_tools_call(request_id: Any, params: Dict[str, Any]) -> None:
    """Handle tools/call request"""
    tool_name = params.get("name")
    arguments = params.get("arguments", {})

    if tool_name not in TOOLS:
        send_response(request_id, error={
            "code": -32601,
            "message": f"Tool not found: {tool_name}"
        })
        return

    try:
        handler = TOOLS[tool_name]["handler"]
        result = handler(arguments)

        # MCP expects content array with text/image objects
        send_response(request_id, {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(result, indent=2)
                }
            ]
        })
    except MCPError as e:
        logger.error(f"Tool error in {tool_name}: {e}")
        send_response(request_id, error={
            "code": -32000,
            "message": str(e)
        })
    except Exception as e:
        logger.exception(f"Unexpected error in {tool_name}")
        send_response(request_id, error={
            "code": -32603,
            "message": f"Internal error: {str(e)}"
        })


def handle_message(msg: Dict[str, Any]) -> None:
    """Route incoming JSON-RPC messages to appropriate handlers"""
    method = msg.get("method")
    request_id = msg.get("id")
    params = msg.get("params", {})

    logger.debug(f"Received: method={method}, id={request_id}")

    if method == "initialize":
        handle_initialize(request_id, params)
    elif method == "tools/list":
        handle_tools_list(request_id)
    elif method == "tools/call":
        handle_tools_call(request_id, params)
    else:
        send_response(request_id, error={
            "code": -32601,
            "message": f"Method not found: {method}"
        })


# ----------------------- MAIN LOOP ----------------------- #

def main() -> None:
    """Main event loop: read JSON-RPC messages from stdin"""
    logger.info(f"Starting {SERVER_INFO['name']} v{SERVER_INFO['version']}")
    logger.info(f"Repository root: {REPO_ROOT}")

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue

        try:
            msg = json.loads(line)
            handle_message(msg)
        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON: {e}")
            # Can't send proper error response without request id
            continue
        except Exception as e:
            logger.exception("Unexpected error processing message")
            continue


if __name__ == "__main__":
    main()

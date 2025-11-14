"""
Simple test client for the nf-treeseq MCP server.

Usage:
    python test_client.py
"""
import json
import subprocess
import sys

def send_request(proc, method, params=None, request_id=1):
    """Send a JSON-RPC request and get response"""
    request = {
        "jsonrpc": "2.0",
        "id": request_id,
        "method": method
    }
    if params:
        request["params"] = params

    request_json = json.dumps(request) + "\n"
    print(f"→ Sending: {request_json.strip()}")
    proc.stdin.write(request_json)
    proc.stdin.flush()

    response_line = proc.stdout.readline()
    print(f"← Received: {response_line.strip()}")
    return json.loads(response_line)

def main():
    # Start the MCP server
    print("Starting MCP server...")
    proc = subprocess.Popen(
        ["python3", "mcp_server/server.py"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    try:
        # 1. Initialize
        print("\n=== Test 1: Initialize ===")
        response = send_request(proc, "initialize", {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {"name": "test-client", "version": "1.0"}
        }, 1)
        assert response["result"]["serverInfo"]["name"] == "nf-treeseq-mcp"
        print("✓ Initialize successful")

        # 2. List tools
        print("\n=== Test 2: List Tools ===")
        response = send_request(proc, "tools/list", request_id=2)
        tools = response["result"]["tools"]
        tool_names = [t["name"] for t in tools]
        print(f"Available tools: {tool_names}")
        assert "listWorkflows" in tool_names
        assert "launchRun" in tool_names
        print("✓ Tools list successful")

        # 3. List workflows
        print("\n=== Test 3: List Workflows ===")
        response = send_request(proc, "tools/call", {
            "name": "listWorkflows",
            "arguments": {}
        }, 3)
        result = json.loads(response["result"]["content"][0]["text"])
        print(f"Found {result['count']} workflows:")
        for wf in result["workflows"][:3]:  # Show first 3
            print(f"  - {wf['name']} ({wf['kind']}): {wf['path']}")
        print("✓ List workflows successful")

        # 4. List runs (should be empty initially)
        print("\n=== Test 4: List Runs ===")
        response = send_request(proc, "tools/call", {
            "name": "listRuns",
            "arguments": {}
        }, 4)
        result = json.loads(response["result"]["content"][0]["text"])
        print(f"Active runs: {result['count']}")
        print("✓ List runs successful")

        # Note: We don't test launchRun here as it would actually start Nextflow
        # To test launchRun manually:
        # response = send_request(proc, "tools/call", {
        #     "name": "launchRun",
        #     "arguments": {
        #         "profile": "test",
        #         "params": {}
        #     }
        # }, 5)

        print("\n=== All tests passed! ===")

    finally:
        proc.terminate()
        proc.wait(timeout=5)

if __name__ == "__main__":
    main()

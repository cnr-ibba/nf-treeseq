#!/usr/bin/env python3
"""
Demo interattivo del server MCP nf-treeseq
Mostra in tempo reale la comunicazione client-server
"""
import subprocess
import json
import time
import sys

def pretty_json(obj):
    """Formatta JSON con colori"""
    return json.dumps(obj, indent=2)

def demo():
    print("=" * 70)
    print("DEMO: Server MCP nf-treeseq")
    print("=" * 70)
    print()

    print("Avvio del server...")
    proc = subprocess.Popen(
        ["python3", "mcp_server/server.py"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    time.sleep(0.5)
    print("✓ Server avviato\n")

    demos = [
        {
            "title": "1. Initialize - Handshake iniziale",
            "request": {
                "jsonrpc": "2.0",
                "id": 1,
                "method": "initialize",
                "params": {
                    "protocolVersion": "2024-11-05",
                    "capabilities": {},
                    "clientInfo": {"name": "demo-client", "version": "1.0"}
                }
            }
        },
        {
            "title": "2. Tools/List - Scopri tool disponibili",
            "request": {
                "jsonrpc": "2.0",
                "id": 2,
                "method": "tools/list"
            }
        },
        {
            "title": "3. ListWorkflows - Trova tutti i workflow",
            "request": {
                "jsonrpc": "2.0",
                "id": 3,
                "method": "tools/call",
                "params": {
                    "name": "listWorkflows",
                    "arguments": {}
                }
            }
        },
        {
            "title": "4. ListRuns - Mostra run attivi",
            "request": {
                "jsonrpc": "2.0",
                "id": 4,
                "method": "tools/call",
                "params": {
                    "name": "listRuns",
                    "arguments": {}
                }
            }
        }
    ]

    try:
        for demo_item in demos:
            print("-" * 70)
            print(f"📤 {demo_item['title']}")
            print("-" * 70)

            request = demo_item['request']
            print("\n🔹 REQUEST:")
            print(pretty_json(request))

            # Invia richiesta
            proc.stdin.write(json.dumps(request) + "\n")
            proc.stdin.flush()

            # Leggi risposta
            response_line = proc.stdout.readline()
            response = json.loads(response_line)

            print("\n🔸 RESPONSE:")

            # Formattazione speciale per alcuni campi
            if "result" in response:
                result = response["result"]

                # Per tools/list, mostra solo nomi
                if "tools" in result and isinstance(result["tools"], list):
                    tool_names = [t.get("name", "?") for t in result["tools"]]
                    print(f'  Tools: {tool_names}')
                    print(f'  (Totale: {len(result["tools"])} tools)')

                # Per listWorkflows, mostra count e primi workflow
                elif "content" in result and len(result["content"]) > 0:
                    text = result["content"][0].get("text", "")
                    try:
                        parsed = json.loads(text)
                        if "workflows" in parsed:
                            print(f'  Workflows trovati: {parsed["count"]}')
                            for wf in parsed["workflows"][:3]:
                                print(f'    - {wf["name"]} ({wf["kind"]})')
                            if parsed["count"] > 3:
                                print(f'    ... e altri {parsed["count"] - 3}')
                        elif "runs" in parsed:
                            print(f'  Run attivi: {parsed["count"]}')
                            if parsed["count"] > 0:
                                for run in parsed["runs"]:
                                    print(f'    - {run["runId"]}: {run["status"]}')
                    except:
                        print(pretty_json(result))
                else:
                    print(pretty_json(result))
            else:
                print(pretty_json(response))

            print()
            time.sleep(1)

        print("=" * 70)
        print("✅ Demo completata!")
        print("=" * 70)
        print()
        print("💡 Per usare il server:")
        print("  1. Test completo: python3 mcp_server/test_client.py")
        print("  2. Avvio manuale: python3 mcp_server/server.py")
        print("  3. Guida VS Code: cat mcp_server/USAGE.md")
        print()

    finally:
        proc.terminate()
        proc.wait(timeout=3)

if __name__ == "__main__":
    demo()

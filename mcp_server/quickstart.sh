#!/bin/bash
# Quick start script per il server MCP nf-treeseq

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "=== nf-treeseq MCP Server - Quick Start ==="
echo ""
echo "Repository: $REPO_ROOT"
echo ""

# Check Python
if ! command -v python3 &> /dev/null; then
    echo "❌ Python3 non trovato. Installa Python 3.8+ prima di continuare."
    exit 1
fi
echo "✓ Python3: $(python3 --version)"

# Check Nextflow
if ! command -v nextflow &> /dev/null; then
    echo "⚠️  Nextflow non trovato nel PATH."
    echo "   Il server funzionerà, ma non potrà lanciare pipeline."
    echo "   Installa Nextflow: https://www.nextflow.io/docs/latest/getstarted.html"
else
    echo "✓ Nextflow: $(nextflow -version 2>&1 | head -n1)"
fi

echo ""
echo "Cosa vuoi fare?"
echo ""
echo "1) Testare il server (client di test)"
echo "2) Avviare il server (modalità interattiva)"
echo "3) Mostrare esempio di configurazione VS Code"
echo "4) Mostrare esempio di chiamata manuale"
echo ""
read -p "Scelta [1-4]: " choice

case $choice in
    1)
        echo ""
        echo "=== Esecuzione test client ==="
        cd "$REPO_ROOT"
        python3 mcp_server/test_client.py
        ;;
    2)
        echo ""
        echo "=== Avvio server (premi Ctrl+C per fermare) ==="
        echo "Invia messaggi JSON-RPC su stdin, ricevi risposte su stdout"
        echo ""
        cd "$REPO_ROOT"
        python3 mcp_server/server.py
        ;;
    3)
        echo ""
        echo "=== Configurazione VS Code ==="
        echo ""
        echo "Aggiungi questo a .vscode/settings.json (o settings utente):"
        echo ""
        cat <<EOF
{
  "mcp.servers": {
    "nf-treeseq": {
      "command": "python3",
      "args": ["$REPO_ROOT/mcp_server/server.py"],
      "env": {},
      "autoStart": true
    }
  }
}
EOF
        echo ""
        echo "Poi ricarica VS Code: Ctrl+Shift+P → 'Developer: Reload Window'"
        ;;
    4)
        echo ""
        echo "=== Esempio chiamata manuale ==="
        echo ""
        echo "# 1. Avvia il server in un terminale:"
        echo "python3 $REPO_ROOT/mcp_server/server.py"
        echo ""
        echo "# 2. In un altro terminale, invia comandi (esempio initialize):"
        echo 'echo '"'"'{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"manual","version":"1.0"}}}'"'"' | python3 '"$REPO_ROOT/mcp_server/server.py"
        echo ""
        echo "# 3. Lista tools disponibili:"
        echo 'echo '"'"'{"jsonrpc":"2.0","id":2,"method":"tools/list"}'"'"' | python3 '"$REPO_ROOT/mcp_server/server.py"
        echo ""
        echo "# 4. Chiama tool listWorkflows:"
        echo 'echo '"'"'{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"listWorkflows","arguments":{}}}'"'"' | python3 '"$REPO_ROOT/mcp_server/server.py"
        ;;
    *)
        echo "Scelta non valida."
        exit 1
        ;;
esac

echo ""
echo "=== Fine ==="

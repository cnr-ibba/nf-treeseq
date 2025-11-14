# Come usare il server MCP nf-treeseq con VS Code

## ✅ Server creato e testato con successo!

Il server MCP è stato creato in `mcp_server/` e include:
- ✅ 4 tools funzionanti (listWorkflows, launchRun, monitorRun, listRuns)
- ✅ Protocollo JSON-RPC 2.0 conforme allo standard MCP
- ✅ Test client funzionante
- ✅ Documentazione completa

## 🚀 Avvio rapido

### Opzione 1: Test immediato
```bash
cd /home/paolo/Projects/NEXTFLOWetude/nf-treeseq
python3 mcp_server/test_client.py
```

### Opzione 2: Script interattivo
```bash
./mcp_server/quickstart.sh
```
Scegli l'opzione 1 per testare o 2 per avviare in modalità interattiva.

## 🔧 Integrazione con VS Code

### Metodo 1: Task VS Code (già configurato!)

Ho creato `.vscode/tasks.json` con due task:

1. **Avviare il server**:
   - Premi `Ctrl+Shift+P`
   - Cerca "Tasks: Run Task"
   - Seleziona "Start nf-treeseq MCP Server"
   - Il server partirà in un terminale dedicato

2. **Testare il server**:
   - Premi `Ctrl+Shift+P`
   - Cerca "Tasks: Run Task"
   - Seleziona "Test MCP Server"

### Metodo 2: Configurazione client MCP (se supportato)

Alcune estensioni VS Code supportano il protocollo MCP direttamente. Per configurarle:

1. Apri le impostazioni VS Code in formato JSON:
   - `Ctrl+Shift+P` → "Preferences: Open User Settings (JSON)"

2. Aggiungi (il nome del campo dipende dall'estensione):

```json
{
  "mcp.servers": {
    "nf-treeseq": {
      "command": "python3",
      "args": ["/home/paolo/Projects/NEXTFLOWetude/nf-treeseq/mcp_server/server.py"],
      "env": {},
      "autoStart": true
    }
  }
}
```

**Oppure**, se usi GitHub Copilot con supporto MCP:
```json
{
  "github.copilot.advanced": {
    "mcpServers": {
      "nf-treeseq": {
        "command": "python3",
        "args": ["/home/paolo/Projects/NEXTFLOWetude/nf-treeseq/mcp_server/server.py"]
      }
    }
  }
}
```

3. Ricarica VS Code:
   - `Ctrl+Shift+P` → "Developer: Reload Window"

### Metodo 3: Estensioni che potrebbero supportare MCP

Cerca nel marketplace VS Code estensioni come:
- "Model Context Protocol" (se disponibile)
- "MCP Client" (se disponibile)
- Alcune versioni di "Claude" o altri AI assistant

**Nota**: Al momento (novembre 2024) il supporto MCP in VS Code è emergente. Se non trovi estensioni native, puoi comunque usare i Task (Metodo 1) per avviare il server e interagire tramite altri client.

## 📱 Uso con altri client MCP

### Claude Desktop (se installato)

1. Crea/modifica `~/.config/claude/config.json`:
```json
{
  "mcpServers": {
    "nf-treeseq": {
      "command": "python3",
      "args": ["/home/paolo/Projects/NEXTFLOWetude/nf-treeseq/mcp_server/server.py"]
    }
  }
}
```

2. Riavvia Claude Desktop

3. Ora puoi chiedere:
   - "Lista i workflow disponibili"
   - "Lancia il workflow di test"
   - "Mostra lo stato dei run"

### Client Python custom

Vedi `mcp_server/test_client.py` come esempio di come creare un client Python.

## 🛠️ Tools disponibili

Una volta connesso, puoi usare questi strumenti in linguaggio naturale:

### 1. **listWorkflows**
"Mostrami tutti i workflow disponibili"
→ Ricevi lista di workflow e subworkflow con path

### 2. **launchRun**
"Lancia il workflow di test con parametro X=Y"
→ Avvia Nextflow e ricevi runId per tracciamento

### 3. **monitorRun**
"Controlla lo stato del run <runId>"
→ Vedi se è in esecuzione, completato o fallito

### 4. **listRuns**
"Mostrami tutti i run attivi"
→ Elenco di tutti i run tracciati

## 📊 Esempio conversazione

**Tu**: "Lista i workflow disponibili"
**AI** (via MCP): → Chiama `listWorkflows`
**Server**: → Risponde con 9 workflow trovati
**AI**: "Ho trovato 9 workflow: treeseq (main), est_sfs, major, reference..."

**Tu**: "Lancia il workflow di test"
**AI** (via MCP): → Chiama `launchRun` con profile="test"
**Server**: → Avvia Nextflow, risponde con runId=xyz
**AI**: "Ho avviato il run xyz. Vuoi monitorarlo?"

**Tu**: "Sì, mostrami lo stato"
**AI** (via MCP): → Chiama `monitorRun` con runId=xyz
**Server**: → Risponde con status, pid, log
**AI**: "Il run è in esecuzione (PID 12345). Ecco gli ultimi log: ..."

## 🔍 Debugging

### Vedere i log del server

Il server logga su stderr. Per catturare i log:

```bash
python3 mcp_server/server.py 2>mcp-server.log
```

Poi in un altro terminale:
```bash
tail -f mcp-server.log
```

### Verificare comunicazione

Puoi inviare messaggi JSON direttamente:

```bash
# Test initialize
echo '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"test","version":"1.0"}}}' | python3 mcp_server/server.py
```

### Problemi comuni

**"nextflow command not found"**
→ Installa Nextflow: `curl -s https://get.nextflow.io | bash`

**Il server non risponde**
→ Controlla che Python 3.8+ sia installato: `python3 --version`

**VS Code non vede il server**
→ Verifica la configurazione in settings.json e ricarica la finestra

## 🚀 Prossimi passi

1. **Prova il test client**: `python3 mcp_server/test_client.py`
2. **Configura il tuo client MCP preferito** (VS Code, Claude, custom)
3. **Estendi il server** con nuovi tools:
   - `runTests`: esegui nf-test
   - `validateSchema`: valida nextflow_schema.json
   - `fetchArtifacts`: raccogli report/trace

Tutta la documentazione è in `mcp_server/README.md`.

## 📚 Risorse

- [Model Context Protocol Docs](https://modelcontextprotocol.io/)
- [JSON-RPC 2.0 Spec](https://www.jsonrpc.org/specification)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)

---

**Il server è pronto all'uso!** 🎉

Inizia con il test client, poi configura il tuo client MCP preferito per sfruttare l'automazione conversazionale delle pipeline Nextflow.

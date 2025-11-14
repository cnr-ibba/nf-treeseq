# nf-treeseq MCP Server

Server MCP (Model Context Protocol) per gestire pipeline Nextflow del progetto nf-treeseq.

## Cosa fa

Questo server espone strumenti (tools) che permettono agli AI assistant di:
- Elencare workflow disponibili
- Lanciare esecuzioni di pipeline Nextflow
- Monitorare lo stato delle esecuzioni
- Visualizzare log e risultati

## Requisiti

- Python 3.8+
- Nextflow installato e disponibile nel PATH
- Repository nf-treeseq

## Installazione

```bash
# Nessuna dipendenza esterna richiesta (solo libreria standard Python)
# Rendi eseguibile il server
chmod +x mcp_server/server.py
```

## Test rapido

```bash
# Test con il client di esempio
cd /home/paolo/Projects/NEXTFLOWetude/nf-treeseq
python3 mcp_server/test_client.py
```

Output atteso:
```
Starting MCP server...
=== Test 1: Initialize ===
✓ Initialize successful
=== Test 2: List Tools ===
Available tools: ['listWorkflows', 'launchRun', 'monitorRun', 'listRuns']
✓ Tools list successful
...
```

## Integrazione con VS Code

### Opzione 1: Usando l'estensione ufficiale MCP (se disponibile)

Se hai installato un'estensione VS Code che supporta MCP (come alcune versioni di GitHub Copilot o Claude):

1. Apri VS Code Settings (JSON): `Ctrl+Shift+P` → "Preferences: Open User Settings (JSON)"

2. Aggiungi la configurazione del server MCP:

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

3. Ricarica VS Code o riavvia la finestra: `Ctrl+Shift+P` → "Developer: Reload Window"

4. Ora puoi usare i tool in conversazione:
   - "Lista i workflow disponibili"
   - "Lancia il workflow di test"
   - "Mostrami lo stato dei run attivi"

### Opzione 2: Avvio manuale per test

Se non hai un'estensione MCP integrata, puoi comunque testare il server manualmente:

1. **Terminale integrato**: Apri un terminale in VS Code e avvia il server:
   ```bash
   python3 mcp_server/server.py
   ```

2. **Invia comandi** (da un altro terminale o usando echo/printf):
   ```bash
   # Initialize
   echo '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"test","version":"1.0"}}}' | python3 mcp_server/server.py

   # List tools
   echo '{"jsonrpc":"2.0","id":2,"method":"tools/list"}' | python3 mcp_server/server.py

   # List workflows
   echo '{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"listWorkflows","arguments":{}}}' | python3 mcp_server/server.py
   ```

### Opzione 3: Task VS Code per avvio rapido

Crea `.vscode/tasks.json`:

```json
{
  "version": "2.0.0",
  "tasks": [
    {
      "label": "Start MCP Server",
      "type": "shell",
      "command": "python3",
      "args": ["${workspaceFolder}/mcp_server/server.py"],
      "isBackground": true,
      "problemMatcher": [],
      "presentation": {
        "reveal": "always",
        "panel": "new"
      }
    }
  ]
}
```

Poi: `Ctrl+Shift+P` → "Tasks: Run Task" → "Start MCP Server"

## Tools disponibili

### 1. listWorkflows

Lista tutti i workflow e subworkflow disponibili.

**Input**: `{}`

**Output**:
```json
{
  "workflows": [
    {
      "name": "nf-treeseq",
      "path": "workflows/nf-treeseq.nf",
      "kind": "workflow"
    }
  ],
  "count": 1
}
```

### 2. launchRun

Lancia un'esecuzione di pipeline Nextflow.

**Input**:
```json
{
  "profile": "test",
  "entry": "workflows/nf-treeseq.nf",
  "params": {
    "input": "assets/samplesheet.csv",
    "outdir": "results"
  }
}
```

**Output**:
```json
{
  "runId": "550e8400-e29b-41d4-a716-446655440000",
  "status": "STARTED",
  "pid": 12345,
  "entry": "workflows/nf-treeseq.nf",
  "profile": "test",
  "startTime": "2024-11-14T10:30:00",
  "runDir": ".runs/550e8400-e29b-41d4-a716-446655440000"
}
```

### 3. monitorRun

Controlla lo stato di un run attivo.

**Input**:
```json
{
  "runId": "550e8400-e29b-41d4-a716-446655440000",
  "includeLogs": true
}
```

**Output**:
```json
{
  "runId": "550e8400-e29b-41d4-a716-446655440000",
  "status": "RUNNING",
  "pid": 12345,
  "startTime": "2024-11-14T10:30:00",
  "workflow": "workflows/nf-treeseq.nf",
  "profile": "test",
  "logs": ["N E X T F L O W  ~  version 23.10.0", "..."]
}
```

### 4. listRuns

Mostra tutti i run tracciati.

**Input**: `{}`

**Output**:
```json
{
  "runs": [
    {
      "runId": "550e8400-e29b-41d4-a716-446655440000",
      "status": "COMPLETED",
      "workflow": "workflows/nf-treeseq.nf",
      "profile": "test",
      "startTime": "2024-11-14T10:30:00",
      "pid": 12345
    }
  ],
  "count": 1
}
```

## Sicurezza

Il server implementa le seguenti misure di sicurezza:

- ✅ **Path confinement**: I path devono essere all'interno della repository
- ✅ **Command whitelist**: Solo `nextflow run` è permesso (no shell arbitraria)
- ✅ **Isolamento run**: Ogni esecuzione ha directory separata in `.runs/`
- ✅ **Process tracking**: I PID dei processi sono tracciati
- ✅ **Logging**: Tutte le operazioni sono loggate su stderr

⚠️ **Per uso in produzione** considera di aggiungere:
- Autenticazione (token/API key)
- Rate limiting
- Timeout per operazioni lunghe
- Gestione memoria/CPU (max run simultanei)
- Audit log persistente

## Struttura file

```
mcp_server/
├── server.py          # Server MCP principale
├── test_client.py     # Client di test
└── README.md          # Questa documentazione

.runs/                 # Directory delle esecuzioni (creata automaticamente)
└── <runId>/
    ├── params.json
    ├── work/
    ├── nextflow.stdout.log
    └── nextflow.stderr.log
```

## Troubleshooting

### "nextflow command not found"
Assicurati che Nextflow sia installato e nel PATH:
```bash
which nextflow
nextflow -version
```

### Il server non risponde
Controlla i log su stderr. Il server logga tutte le operazioni:
```bash
python3 mcp_server/server.py 2>debug.log
```

### Run non parte
Verifica i parametri e che il workflow esista:
```bash
ls workflows/nf-treeseq.nf
```

### VS Code non vede il server
1. Controlla che la configurazione sia in `settings.json`
2. Verifica che il path assoluto sia corretto
3. Ricarica la finestra VS Code
4. Guarda Output → "MCP" per eventuali errori

## Prossimi sviluppi

- [ ] Tool `validateSchema`: validazione nextflow_schema.json
- [ ] Tool `runTests`: esecuzione suite nf-test
- [ ] Tool `fetchArtifacts`: raccolta report/trace/timeline
- [ ] Streaming log in tempo reale (SSE)
- [ ] Dashboard web per monitoraggio
- [ ] Integrazione con Tower/Seqera Platform

## Licenza

Stesso della repository nf-treeseq (vedi LICENSE principale).

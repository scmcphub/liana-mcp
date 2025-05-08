# Liana-MCP

Natural language interface for scRNA-Seq analysis with Liana through MCP.

## 🪩 What can it do?

- IO module like read and write scRNA-Seq data
- Cell-cell communication method
- Plotting module, circle plot, dotplot

## ❓ Who is this for?

- Anyone who wants to do scRNA-Seq analysis natural language!
- Agent developers who want to call Liana's functions for their applications

## 🌐 Where to use it?

You can use Liana-mcp in most AI clients, plugins, or agent frameworks that support the MCP:

- AI clients, like Cherry Studio
- Plugins, like Cline
- Agent frameworks, like Agno 

## 🎬 Demo

A demo showing scRNA-Seq cell cluster analysis in a AI client Cherry Studio using natural language based on Liana-mcp


## 🏎️ Quickstart

### Install

Install from PyPI
```
pip install Liana-mcp
```
you can test it by running
```
liana-mcp run
```

#### run scnapy-server locally
Refer to the following configuration in your MCP client:

```
"mcpServers": {
  "liana-mcp": {
    "command": "liana-mcp",
    "args": [
      "run"
    ]
  }
}
```

#### run scnapy-server remotely
Refer to the following configuration in your MCP client:

run it in your server
```
Liana-mcp run --transport shttp --port 8000
```

Then configure your MCP client, like this:
```
http://localhost:8000/mcp
```

## 🤝 Contributing

If you have any questions, welcome to submit an issue, or contact me(hsh-me@outlook.com). Contributions to the code are also welcome!

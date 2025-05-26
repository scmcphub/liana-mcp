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

https://github.com/user-attachments/assets/40fb5bd8-a166-4993-9979-3258ef6646a0


## 📚 Documentation

scmcphub's complete documentation is available at https://docs.scmcphub.org

## 🏎️ Quickstart

### Install

Install from PyPI
```
pip install liana-mcp
```
you can test it by running
```
liana-mcp run
```


#### run liana-mcp locally
Refer to the following configuration in your MCP client:

check path
```
$ which liana 
/home/test/bin/liana-mcp
```

```
"mcpServers": {
  "liana-mcp": {
    "command": "/home/test/bin/liana-mcp",
    "args": [
      "run"
    ]
  }
}
```

#### run liana-server remotely
Refer to the following configuration in your MCP client:

run it in your server
```
liana-mcp run --transport shttp --port 8000
```

Then configure your MCP client in local AI client, like this:
```

"mcpServers": {
  "liana-mcp": {
    "url": "http://localhost:8000/mcp"
  }
}
```

## 🤝 Contributing

If you have any questions, welcome to submit an issue, or contact me(hsh-me@outlook.com). Contributions to the code are also welcome!

## Citing

If you use liana-mcp in for your research, please consider citing  following work: 
> Dimitrov D., Schäfer P.S.L, Farr E., Rodriguez Mier P., Lobentanzer S., Badia-i-Mompel P., Dugourd A., Tanevski J., Ramirez Flores R.O. and Saez-Rodriguez J. LIANA+ provides an all-in-one framework for cell–cell communication inference. Nat Cell Biol (2024). https://doi.org/10.1038/s41556-024-01469-w

"""
Command-line interface for liana-mcp.

This module provides a CLI entry point for the liana-mcp package.
"""
from scmcp_shared.cli import MCPCLI
from .server import LianaMCPManager

cli = MCPCLI(
    name="liana-mcp", 
    help_text="Liana MCP Server CLI",
    manager=LianaMCPManager
)


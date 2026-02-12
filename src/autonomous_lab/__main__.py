#!/usr/bin/env python3
"""
Autonomous Lab - Main entry point
"""

import argparse
import sys


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Autonomous Lab - MCP server for iterative research sessions"
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    subparsers.add_parser("server", help="Start MCP server (default)")
    subparsers.add_parser("version", help="Show version")

    args = parser.parse_args()

    if args.command == "version":
        from . import __author__, __version__
        print(f"Autonomous Lab v{__version__}")
        print(f"Author: {__author__}")
    elif args.command == "server" or args.command is None:
        from .server import main as server_main
        return server_main()
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()

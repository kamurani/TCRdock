import click as ck 

from .setup import commands as setup_commands

@ck.group()
def entry_point():
    pass

entry_point.add_command(setup_commands.setup)

if __name__ == "__main__":
    entry_point()
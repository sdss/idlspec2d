from __future__ import annotations

from build_doc_legacy import build_legacy_docs

from pathlib import Path
import click
import subprocess
import shutil

try:
    from pkg_resources import resource_filename
    bindir = Path(resource_filename("boss_drp", "../../bin/"))
    prodir = Path(resource_filename("boss_drp", "../../pro/"))
    docdir = Path(resource_filename("boss_drp", "../../docs/sphinx/"))
except Exception:
    file_path = Path(__file__).resolve()
    bindir = file_path.parent.parent / "bin"
    prodir = file_path.parent.parent / "pro"
    docdir = file_path.parent.parent / "docs" / "sphinx"

HEADING_CHARS = ["=", "-", "^", '"', "~", "+", "#", "*"]


def _safe_ctx(cmd: click.Command, info_name: str) -> click.Context:
    return click.Context(cmd, info_name=info_name)


def rst_heading(title: str, level: int, toc: bool = True) -> str:
    ch = HEADING_CHARS[min(level, len(HEADING_CHARS) - 1)]
    parts = []
    parts.append(f"{title}\n{ch * len(title)}\n\n")     
    if toc: 
        parts.append('.. contents::\n')
        parts.append('    :depth: 3'+'\n')
        parts.append('    :local:\n')
        parts.append('    :class: this-will-duplicate-information-and-it-is-still-useful-here\n')
        parts.append('    :backlinks: none\n\n')
    return ''.join(parts)


def list_subcommands(cmd: click.Command) -> list[tuple[str, click.Command]]:
    if not isinstance(cmd, click.Group):
        return []

    ctx = _safe_ctx(cmd, info_name=cmd.name or "cli")
    try:
        items: list[tuple[str, click.Command]] = []
        for sub_name in cmd.list_commands(ctx):
            sub_cmd = cmd.get_command(ctx, sub_name)
            if sub_cmd is not None:
                items.append((sub_name, sub_cmd))
        return items
    finally:
        ctx.close()


def has_subcommands(cmd: click.Command) -> bool:
    return bool(list_subcommands(cmd))


def render_help_dropdown(name: str, cmd: click.Command, title: str,
                         collapse: bool = True, main: bool = False) -> str:
    ctx = _safe_ctx(cmd, info_name=name)#.split()[-1])
    try:
        help_text = cmd.get_help(ctx).rstrip()
    finally:
        ctx.close()
    
    if main:
        return (
            f".. _{'_'.join(name.split())}_py:\n\n"
            f".. code-block:: text\n\n"
            + "\n".join(f"   {line}" for line in help_text.splitlines())
            + "\n\n"
        )

    return (
        f".. _{'_'.join(name.split())}_py:\n\n"
        f".. admonition:: {title}\n"
        f"   :collapsible: {'closed' if collapse else 'open'}\n\n"
        f"   .. code-block:: text\n\n"
        + "\n".join(f"      {line}" for line in help_text.splitlines())
        + "\n\n"
    )


def render_command_tree(cmd: click.Command, name: str, depth: int = 0, path: tuple[str, ...] = (), func: str = 'cli') -> str:
    parts: list[str] = []

    current_path = path + ((cmd.name or name),)
    full_name = " ".join(current_path).strip()

    
    title = list(current_path)
    title[0] = name
    title = ' '.join(title).strip()

    if depth == 0:
        parts.append(f'.. _{name}:\n\n')
        parts.append(rst_heading(name, 2))
        parts.append(render_help_dropdown(name, cmd, title, main=True))
    elif has_subcommands(cmd):
        parts.append(rst_heading(title, depth+2))
        parts.append(render_help_dropdown(full_name, cmd, title, collapse= False))
    else:
        parts.append(render_help_dropdown(full_name, cmd, title))

    for sub_name, sub_cmd in list_subcommands(cmd):
        if current_path[-1] == func:
             current_path =  current_path[:-1]
        
        parts.append(render_command_tree(sub_cmd, name, depth + 1, func = func, path = current_path))# + (sub_name,)))

    return "".join(parts)


def export_click_help_to_rst(cli: click.Command, name: str, func: str = 'cli') -> list[str]:
    print(f"Building Docs for {name}")
    parts: list[str] = []
    parts.append(render_command_tree(cli, name, path = (name,), func = func))
    return parts



def build_docs():
    # Import your CLIs here
    from boss_drp.cli.BOSS_drp import cli as boss_drp_cli
    from boss_drp.cli.SOS import cli as sos_cli
    from boss_drp.cli.flatlib import cli as flatlib_cli
    from boss_drp.cli.boss_drp.version import version as version_cli
    from boss_drp.cli.boss_drp.tools import run_sdR_hdrfix
    from boss_drp.cli.SOS import run_log


    parts: list[str] = []
    parts.append(":tocdepth: 5\n\n")
    parts.append(".. highlight:: none\n\n")
    parts.append(rst_heading('Full Command Documention', level = 0, toc= False))
    parts.append('Documented below are the primary commands used to run the BOSS Data Reduction Pipeline. '+
                 'However, there are numerous other routines included in this package, '+
                 'which are called by these commands and have their own internal documentation.'+
                 'The legacy CLI interface is still included, documneted on :doc:`Legacy CLI<doc_legacy>`\n\n')

    parts.append(rst_heading('Full Python Command Usage', 1))
    parts.extend(export_click_help_to_rst(boss_drp_cli, name="boss_drp", func = 'cli'))
    parts.extend(export_click_help_to_rst(sos_cli, name="SOS"))
    parts.extend(export_click_help_to_rst(flatlib_cli, name="boss_flatlib"))
    parts.extend(export_click_help_to_rst(version_cli, name="idlspec2d_version"))
    parts.extend(export_click_help_to_rst(run_sdR_hdrfix, name="sdR_hdrfix"))
    parts.extend(export_click_help_to_rst(run_log, name="BOSS_log"))

    parts.append(rst_heading('Full Bash Command Usage', 1))
    for command in sorted(bindir.glob("*.bash"), key=lambda p: p.name):
        print(f"Building Docs for {command.name}")
        docstr = subprocess.getoutput(f'{command} -h')
        parts.append(rst_heading(command.name,2))
        docstr = "".join("   " + line for line in docstr.splitlines(True))
        parts.append('::\n\n')
        parts.extend(docstr+'\n\n')

    
    parts.append(rst_heading('IDL Command Usage', 1))
    for command in ['spreduce2d.pro','rm_combine_script.pro',
                    'spreduce1d_empca.pro',
                    'spspec_target_merge.pro']:
        print(f"Building Docs for {command}")

        pf = list(prodir.glob(f"*/{command}"))
        if len(pf) == 0: continue

        docstr = []
        with open(pf[0],'r') as prof:
            docstr = prof.read()
        dss = []
        for ds in docstr.split('\n'):
            if ';------------------' in ds:
                break
            if len(ds.strip()) == 0:
                continue
            dss.append("   " +ds)
        docstr = '\n'.join(dss)
        parts.append(rst_heading(command,2))
        parts.append('::\n\n')
        parts.extend(docstr+'\n\n')

    parts.append('.. highlight:: defaults\n')
    parts.append('\n.. End of document\n')
    out_file = docdir / "doc.rst"
    out_file.write_text("".join(parts), encoding="utf-8")
    print(f"Wrote {out_file}")


    # Copy Tree schema diagrams for readthedocs
    src_dir = docdir.resolve().parent.parent / "datamodel" / "tree"
    dst_dir = docdir / "_static" / "tree"

    dst_dir.mkdir(parents=True, exist_ok=True)

    print(src_dir)
    print(list(src_dir.glob("*.png")))

    for img in src_dir.glob("*.png"):
        print(f"Copying {img} to {dst_dir}")
        shutil.copy(img, dst_dir)

if __name__ == '__main__' :
    """
    Build BOSS DRP Documention
    """
    build_docs()
    build_legacy_docs()

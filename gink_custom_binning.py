import os, re

TOOL = 'gink_custom_binning'
PREFIX = TOOL + '_'
DEFAULT = 'variable_175000_48_bwa'
VALID = re.compile(r'^[A-Za-z0-9_.-]+$')

def setup(d4, tools, binnings):
    tools, binnings = list(tools), list(dict.fromkeys(binnings))
    if TOOL in tools and not binnings: raise ValueError(F'{TOOL} requires --binnings')
    tools = [tool for tool in tools if tool != TOOL]
    if not binnings: return tools
    for binning in binnings:
        if not VALID.fullmatch(binning): raise ValueError(F'Invalid Ginkgo binning: {binning}')
        tool = PREFIX + binning
        if tool not in tools: tools.append(tool)
        d4.SC_CN_TOOL_DEPENDENCY_TO_DEPENDENT['bam2bed'][tool] = ''
        d4.SC_CN_TOOL_TO_RUN_ORDER[tool] = d4.SC_CN_TOOL_TO_RUN_ORDER['ginkgo']
        d4.SC_CN_TOOL_TO_RUN_MODE[tool] = d4.SC_CN_TOOL_TO_RUN_MODE['ginkgo']
        d4.SC_CN_EVAL_TOOLS.add(tool)
    if 'bam2bed' not in tools: tools.append('bam2bed')
    if not getattr(d4.run_tool_1, '_gink_custom_binning', False): d4.run_tool_1 = _wrap(d4.run_tool_1)
    return tools

def _wrap(original):
    def run_tool_1(infodict, tool, *args, **kwargs):
        if not tool.startswith(PREFIX): return original(infodict, tool, *args, **kwargs)
        binning = tool[len(PREFIX):]
        ret = original(infodict, 'ginkgo', *args, **kwargs)
        script = kwargs.get('script', args[2] if len(args) > 2 else None)
        if script and os.path.exists(script):
            text = open(script).read().replace(F'--binning {DEFAULT}', F'--binning {binning}').replace('#sequential=run.ginkgo/', F'#sequential=run.{tool}/')
            open(script, 'w').write(text)
        deps, cmds, bam2bed, lib2bed = ret
        deps = [tuple((x.replace('_tool_ginkgo.rule', F'_tool_{tool}.rule') if isinstance(x, str) else x) for x in dep) for dep in deps]
        cmds = [cmd.replace(F'--binning {DEFAULT}', F'--binning {binning}').replace('#sequential=run.ginkgo/', F'#sequential=run.{tool}/') for cmd in cmds]
        return deps, cmds, bam2bed, lib2bed
    run_tool_1._gink_custom_binning = True
    return run_tool_1

from constants import REPORTER_LINES, DRIVER_LINES, FLUOROPHORES


def parse_genotype(gtype):
    """Return a dict of information about the drivers and reporters in a
    genotype string like:
    
        'Tlx3-Cre_PL56/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt'
        
    Returns
    -------
    drivers : list
        Transgenic drivers in this genotype (eg: pvalb, sst, etc.)
    reporters : list
        Transgenic reporters in this genotype (eg: tdTomato, EGFP, etc.)
    reporter_map : dict
        Dictionary that maps from driver to reporter (forward and backward)
    color_map : dict
        Dictionary that maps from driver to reporter color (forward and backward)
    driver_lines : list
        Transgenic driver lines in this genotype
    reporter_lines : list
        Transgenic reporter lines in this genotype
    """
    ignore = ['wt', 'PhiC31-neo']
    
    parts = set()
    for part in gtype.split(';'):
        for subpart in part.split('/'):
            if subpart in ignore:
                continue
            parts.add(subpart)
            
    driver_lines = [p for p in parts if p in DRIVER_LINES]
    reporter_lines = [p for p in parts if p in REPORTER_LINES]
    extra = parts - set(driver_lines+reporter_lines)
    if len(extra) > 0:
        raise Exception("Unknown genotype part(s): %s" % str(extra))
    
    # map {cre line : fluorophore : color} and back again
    reporter_map = {}
    color_map = {}
    drivers = set()
    reporters = set()
    
    for d in driver_lines:
        d_creflp, driver = DRIVER_LINES[d]
        drivers.add(driver)
        for r in reporter_lines:
            r_creflp, reporter = REPORTER_LINES[r]
            if r_creflp != d_creflp:
                # this is an oversimplification, but should work for now..
                continue
            
            reporter_map[driver] = reporter
            reporter_map[reporter] = driver
            reporters.add(reporter)
            
            color = FLUOROPHORES[reporter]
            color_map[driver] = color
            color_map[color] = driver
    
    return {'driver_lines': driver_lines, 'reporter_lines': reporter_lines,
            'drivers': list(drivers), 'reporters': list(reporters), 
            'reporter_map': reporter_map, 'color_map': color_map}
    
    

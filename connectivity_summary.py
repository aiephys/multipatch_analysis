"""
I have a summary of ALM multipatch experiments stored on workflowy and exported to expt_summary.txt.

This file describes, for each experiment:
   - which cre labels were used
   - which cells are +/- for which labels
   - which cells had usable recordings
   - which cells were connected

The purpose of this script is to parse the exported file and return a summary
of total connections detected / probed for each cell type pair.

"""
import re, sys, traceback


ALL_CRE_TYPES = ['sst', 'pvalb', 'tlx3', 'sim1', 'rorb']

# read exported file

#steph = open('/home/luke/mnt/mp1/Steph/DataSummary', 'r').readlines()
#pasha = open('/home/luke/mnt/mp2/data/Pasha/connection_analysis', 'r').readlines()
#alex = open('/home/luke/mnt/mp3/data/Alex/connection analysis', 'r').readlines()
#lines = steph + pasha + alex
lines = open(sys.argv[1], 'r').readlines()


# first parse indentation to generate a hierarchy of strings
# Note: summary was exported in plain text mode, which can be difficult to parse.
# Might want to switch to XML if this becomes an issue.

class Entry(object):
    def __init__(self, line, parent):
        self.indentation = indentation(line)
        self.lines = []
        self.add_line(line)
        self.parent = parent
        self.children = []
        if parent is not None:
            parent.add_child(self)
        
    def add_line(self, line):
        if indentation(line) != self.indentation:
            raise IndentationError(line)
        self.lines.append(line[self.indentation:].rstrip())
        
    def add_child(self, child):
        self.children.append(child)

    def print_tree(self):
        print("\n".join([('    '*self.indentation)+l for l in self.lines]))
        for ch in self.children:
            ch.print_tree()

        
def indentation(line):
    return len(line) - len(line.lstrip('- '))


root = Entry('', None)
root.indentation = -1
current = root

for line in lines:
    if line.lstrip().startswith('#'):
        continue
    
    if line.strip() == '':
        continue
    ind = indentation(line)
    if ind > current.indentation:
        ch = Entry(line, parent=current)
        current = ch
    else:
        while current.indentation > ind:
            current = current.parent
        ch = Entry(line, parent=current.parent)
        current = ch
        continue

        

# Now go through each experiment and read cell type / connectivity data

class Experiment(object):
    def __init__(self, entry):
        try:
            self.entry = entry
            self.expt_id = entry.lines[0]
            self.cells = {x:Cell(self,x) for x in range(1,9)}
            self.connections = []
            self._summary = None
        except Exception as exc:
            raise Exception("Error parsing experiment: %s" % entry.lines) from exc
        
        for ch in entry.children:
            try:
                if ch.lines[0] == 'Labeling':
                    self.parse_labeling(ch)
                elif ch.lines[0] == 'Cell QC':
                    self.parse_qc(ch)
                elif ch.lines[0] == 'Connections':
                    self.parse_connections(ch)
                elif ch.lines[0] == 'Conditions':
                    continue
                else:
                    raise Exception("Invalid experiment entry %s" % ch.lines[0])
            except Exception as exc:
                raise Exception("Error parsing %s for experiment %s." % (ch.lines[0], self.expt_id)) from exc
            
        # gather a list of all cre types
        cre_types = set()
        for cell in self.cells.values():
            cre_types |= set(cell.labels.keys())
        cre_types.discard('biocytin')
        self.cre_types = list(cre_types)
                
    def parse_labeling(self, entry):
        """
        "Labeling" section should look like:
        
            Labeling:
                sim1: 1+ 2- 3x 4x+ 5+? 6?
                biocytin: ...
                
        This example has the following interpretation:
        
            1+   Cell 1 is reporter-positive
            2-   Cell 2 is reporter-negative    
            3x   Cell 3 type cannot be determined (no pipette tip found)
            4x+  Cell 4 was not dye filled, but pipette tip appears to be touching cre-positive cell
            5+?  Cell 5 looks like cre-positive, but image is somewhat unclear
            6?   Cell 6 is dye filled, but cre type is ambiguous
        """
        for ch in entry.children:
            line = ch.lines[0]
            # line looks like "sim1: 1+ 2-'
            
            parts = re.split('\s+', line)
            
            # first part is cre type and a colon
            assert parts[0].endswith(':')
            cre = parts[0][:-1]
            if not (cre == 'biocytin' or cre in ALL_CRE_TYPES):
                raise Exception("Invalid cre type: %s" % cre)
            
            # parse the remainder of the line
            if len(parts[1:]) == 1 and parts[1].strip() == '?':
                # no data
                continue

            for part in parts[1:]:
                m = re.match('(\d+)(x)?(\+|\-)?(\?)?', part)
                if m is None:
                    raise Exception('invalid label record: %s' % part)
                grps = m.groups()
                cell_id = int(grps[0])
                cell = self.cells[cell_id]
                absent = grps[1] == 'x'
                positive = grps[2]
                uncertain = grps[3] == '?'
                
                assert cre not in cell.labels
                cell.labels[cre] = positive

    def parse_qc(self, entry):
        """Parse cell quality control. Looks like:
        
            Holding: 1- 2+ 3- 4- 6/ 7+
            Access: 1- 2/ 3+ 4- 6/ 7/
            Spiking: 1- 2+ 3+ 4- 6+ 7+
        
        Where + means pass, / means borderline pass, - means fail, ? means unknown
        """
        for ch in entry.children:
            parts = re.split('\s+', ch.lines[0].strip())
            for part in parts[1:]:
                m = re.match(r'(\d+)(\+|\/|\-|\?)', part)
                if m is None:
                    raise Exception('Invalid cell QC string "%s"' % part)
                cell_id = int(m.groups()[0])
                val = m.groups()[1]
                cell = self.cells[cell_id]
                
                if parts[0] == 'Holding:':
                    cell.holding_qc = val in '+/'
                elif parts[0] == 'Access:':
                    cell.access_qc = val in '+/'                    
                elif parts[0] == 'Spiking:':
                    cell.spiking_qc = val in '+/'                    
                else:
                    raise Exception("Invalid Cell QC line: %s" % ch.lines[0])
    
    def parse_connections(self, entry):
        if len(entry.children) == 0 or entry.children[0].lines[0] == 'None':
            return
        for con in entry.children:
            m = re.match(r'(\d+)\s*->\s*(\d+)(.*)', con.lines[0].strip())
            if m is None:
                raise Exception("Invalid connection: %s" % con.lines[0])
            
            if m.groups()[2].strip() == '?':
                # ignore questionable connections
                continue
            self.connections.append((int(m.groups()[0]), int(m.groups()[1])))

    def summary(self):
        """Return a structure summarizing (non)connectivity in the experiment.
        
        Looks like:
        
            {(pre_type, post_type): [connected, unconnected], ...}
        """
        if self._summary is None:
            csum = {}
            for i,ci in self.cells.items():
                for j,cj in self.cells.items():
                    if i == j:
                        continue
                    if ci.spiking_qc is not True:
                        # presynaptic cell failed spike QC; ignore
                        continue
                    if cj.pass_qc is not True:
                        # postsynaptic cell failed ephys qc; ignore
                        continue
                    if ci.cre_type is None or cj.cre_type is None:
                        # indeterminate cell types; ignore
                        #print("Ignore probe (ind. cell type) %s:%d-%d" % (self.expt_id, i, j))
                        #if (i, j) in self.connections:
                            #print("    --> connected!")
                        continue
                    typ = (ci.cre_type, cj.cre_type)
                    if typ not in csum:
                        csum[typ] = [0, 0]
                    if (i, j) in self.connections:
                        csum[typ][0] += 1
                    else:
                        csum[typ][1] += 1
            self._summary = csum
        return self._summary
    
    @property
    def n_connections_probed(self):
        return sum([x[0]+x[1] for x in self.summary().values()])

    @property
    def n_connections(self):
        return sum([x[0] for x in self.summary().values()])
        
        
class Cell(object):
    def __init__(self, expt, cell_id):
        self.expt = expt
        self.cell_id = cell_id
        self.access_qc = None
        self.holding_qc = None
        self.spiking_qc = None
        self.labels = {}
    
    @property
    def pass_qc(self):
        """True if cell passes QC.
        """
        if self.access_qc is True and self.holding_qc is True:
            return True
        elif self.access_qc is False or self.holding_qc is False:
            return False
        
        # None means cell is not present in ephys data
        return None

    @property
    def cre_type(self):
        """Cre type string for this cell.
        
        If the cell is reporter-negative then cre_type is 'unk'.
        If the cell has ambiguous or missing data then cre_type is None.
        """
        default = 'unknown'
        ct = None
        for label,pos in self.labels.items():
            if label == 'biocytin':
                continue
            if pos == '+':
                if ct not in (None, default):
                    raise Exception("%s has multiple labels!" % self)
                ct = label
            elif pos == '-':
                if ct is not None:
                    continue
                ct = default
        return ct

    def __repr__(self):
        return "<Cell %s:%d>" % (self.expt.expt_id, self.cell_id)

        
# Parse experiment data
expts = []
errs = []
for entry in root.children:
    try:
        expts.append(Experiment(entry))
    except Exception as exc:
        errs.append((entry, sys.exc_info()))

if len(errs) > 0:
    print("Errors loading %d experiments:" % len(errs))
    for entry, exc in errs:
        print("=======================")
        print("\n".join(entry.lines))
        traceback.print_exception(*exc)
        print("")

expt_ids = {e.expt_id:e for e in expts}

# Generate summary of experiments
print("----------------------------------------------------------")
print("  Experiment Summary  (# probed, # connected, cre types)")
print("----------------------------------------------------------")

tot_probed = 0
tot_connected = 0
for expt in expts:
    n_p = expt.n_connections_probed
    n_c = expt.n_connections
    tot_probed += n_p
    tot_connected += n_c
    print("%s:  \t%d\t%d\t%s" % (expt.expt_id, n_p, n_c, ', '.join(expt.cre_types)))
print("")


# Generate a summary of connectivity
print("-------------------------------------------------------------")
print("     Connectivity  (# connected/probed, % connectivity)")
print("-------------------------------------------------------------")
summary = {}
for expt in expts:
    for k,v in expt.summary().items():
        if k not in summary:
            summary[k] = [0, 0]
        summary[k][0] += v[0]
        summary[k][1] += v[1]

totals = []
for k,v in summary.items():
    totals.append((k[0], k[1], v[0], v[0]+v[1], 100*v[0]/(v[0]+v[1])))
    
colsize = max([len(t[0]) + len(t[1]) for t in totals]) + 8
totals.sort(key=lambda x: (x[4], x[3]), reverse=True)
for tot in totals:
    pad = " " * (colsize - (len(tot[0]) + len(tot[1]) + 3))
    fields = list(tot)
    fields.insert(2, pad)
    fields = tuple(fields)
    try:
        print("%s → %s%s\t:\t%d/%d\t%0.2f%%" % fields)
    except UnicodeEncodeError:
        print("%s - %s%s\t:\t%d/%d\t%0.2f%%" % fields)

print("\nTotal:  \t%d/%d\t%0.2f%%" % (tot_connected, tot_connected+tot_probed, 100*tot_connected/(tot_connected+tot_probed)))
print("")



print("-----------------------")
print("       Labeling")
print("-----------------------")

n_qc_passed = 0
n_dye_passed = 0
n_qc_and_biocytin_passed = 0
n_dye_and_biocytin_passed = 0

for expt in expts:
    for cell in expt.cells.values():
        biocytin = cell.labels.get('biocytin', None)
        if biocytin is None:
            # ignore cells with no biocytin data
            continue
        biocytin = (biocytin == '+')
        if cell.pass_qc:
            n_qc_passed += 1
            if biocytin:
                n_qc_and_biocytin_passed += 1
        if cell.cre_type is not None:
            n_dye_passed += 1
            if biocytin:
                n_dye_and_biocytin_passed += 1

dye_biocytin_percent = 100*n_dye_and_biocytin_passed/n_dye_passed if n_dye_passed > 0 else 0
print("%0.2f (%d/%d) of dye-filled cells had a biocytin fill" % (dye_biocytin_percent, n_dye_and_biocytin_passed, n_dye_passed))

qc_biocytin_percent = 100*n_qc_and_biocytin_passed/n_qc_passed if n_qc_passed > 0 else 0
print("%0.2f (%d/%d) of qc-passed cells had a biocytin fill" % (qc_biocytin_percent, n_qc_and_biocytin_passed, n_qc_passed))

print("")







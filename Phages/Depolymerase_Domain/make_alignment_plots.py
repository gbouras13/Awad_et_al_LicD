


import pytrimal

import rich.console
import rich.panel
from rich_msa import RichAlignment



alignment = "domains.msa"

ali = pytrimal.Alignment.load(alignment)

# trim to get only the
# order so 10, 16, 20
ali_trimmed = pytrimal.Alignment(
    names=[ali.names[0], ali.names[1], ali.names[2]],
    sequences=[
        ali.sequences[0][684:804],
        ali.sequences[1][684:804],
        ali.sequences[2][684:804]
        ]
        )

# trim to get from 




def show_alignment(alignment, output_name):
    
    console = rich.console.Console(width=len(alignment.sequences[0])+40,record=True)
    widget = RichAlignment(names=[n.decode() for n in alignment.names], sequences=alignment.sequences, max_name_width=30)
    panel = rich.panel.Panel(widget, title_align="left", title="({} residues, {} sequences)".format(len(alignment.sequences[0]), len(alignment.sequences)))
    console.print(panel)


    with open(output_name, "w") as f:
        f.write(console.export_svg())


show_alignment(ali, "full_alignment.svg")

show_alignment(ali_trimmed, "depolymerase_alignment.svg")
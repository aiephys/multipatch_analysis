from multipatch_analysis.database.database import declarative_base, make_table
from multipatch_analysis.database import Database
from multipatch_analysis import config

DataNotesORMBase = declarative_base()


PairNotes = make_table(
    name='pair_notes',
    comment= "User-curated ",
    columns=[
        ('expt_id', 'str', 'Unique experiment identifier (acq_timestamp)', {'index': True}),
        ('pre_cell_id', 'str', 'external id of presynaptic cell'),
        ('post_cell_id', 'str', 'external id of postsynaptic cell'),
        ('notes','object', 'pair data dict which includes synapse call, initial fit parameters, output fit parameters, comments, etc'),
    ],
    ormbase=DataNotesORMBase,
)


db = Database(config.synphys_db_host, config.synphys_db_host_rw, "data_notes", DataNotesORMBase)

db.create_tables()





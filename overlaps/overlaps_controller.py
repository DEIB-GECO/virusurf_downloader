import sys
from overlaps import multi_database_manager


def run():
    inp = sys.argv[7].lower()
    to_commit_str = sys.argv[8].lower()
    if to_commit_str == 'true' or to_commit_str == 'yes':
        multi_database_manager.user_asked_to_commit = True
    else:
        multi_database_manager.user_asked_to_commit = False


    if inp == 'genbank_gisaid':
        from overlaps.genbank_gisaid.genbank_gisaid import run
        run()
    elif inp == 'genbank_nmdc':
        from overlaps.genbank_nmdc.genbank_nmdc import run
        run()
    elif inp == 'gisaid_nmdc':
        from overlaps.gisaid_nmdc.gisaid_nmdc import run
        run()
    elif inp == 'coguk_nmdc':
        from overlaps.coguk_nmdc.coguk_nmdc import run
        run()
    elif inp == 'coguk_gisaid':
        from overlaps.coguk_gisaid.coguk_gisaid import run
        run()
    elif inp == 'coguk_genbank':
        from overlaps.coguk_genbank.coguk_genbank import run
        run()
    elif inp == 'all':
        from overlaps.genbank_gisaid import genbank_gisaid
        genbank_gisaid.run()
        from overlaps.genbank_nmdc import genbank_nmdc
        genbank_nmdc.run()
        from overlaps.gisaid_nmdc import gisaid_nmdc
        gisaid_nmdc.run()
        from overlaps.coguk_nmdc import coguk_nmdc
        coguk_nmdc.run()
        from overlaps.coguk_gisaid import coguk_gisaid
        coguk_gisaid.run()
        from overlaps.coguk_genbank import coguk_genbank
        coguk_genbank.run()
    else:
        raise ValueError('Please specify which data source you wish to analyse for checking overlapping sequences.\n'
                         'Valid choices are:\n'
                         'genbank_gisaid\n'
                         'genbank_nmdc\n'
                         'gisaid_nmdc\n'
                         'coguk_gisaid\n'
                         'coguk_genbank\n'
                         'all (all pairs of data sources will be analyzed, one pair at a time).')
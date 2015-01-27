from distutils.core import setup

setup(
	name='BioClasses',
	version='0.1.5',
	packages=["BioClasses"],
	data_files=[("genetic_codes", \
["data/genetic_codes/euplotid_genetic_code.txt", \
"data/genetic_codes/human_genetic_code.txt", \
"data/genetic_codes/tetrahymena_genetic_code.txt"]), \
("CAI_tables", \
["data/CAI_tables/homo_CAI_table.txt", \
"data/CAI_tables/euplotid_CAI_table.txt", \
"data/CAI_tables/tetrahymena_CAI_table.txt"]), \
("transition_matrices", \
["data/transition_matrices/homo_transition_matrix.pic", \
"data/transition_matrices/euplotid_transition_matrix.pic", \
"data/transition_matrices/tetrahymena_transition_matrix.pic"])],
#	package_data={"BioClasses": ["data/genetic_codes/*.txt", "data/CAI_tables/*.txt", "data/transition_matrices/*.pic"]},	
	requires=["pysam", "Biopython", "scipy"],
	author="Paul K. Korir",
	author_email="paul.korir@gmail.com",
	url="http://www.paulkorir.com/projects/BioClasses",
	license="LICENSE.txt",
	description="Python classes for bioinformatics",
	long_description=open( "README.txt" ).read(),
)

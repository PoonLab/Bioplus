from setuptools import setup

setup(
        name='Bioplus',
        description='in-house scripts that extend the functionality of Biopython',
        py_modules=[
            'Bioplus.add_dates', 
            'Bioplus.codon_align',
            'Bioplus.consensus',
            'Bioplus.dist',
            'Bioplus.filter_fasta',
            'Bioplus.get_metadata',
            'Bioplus.iqtree_asr',
            'Bioplus.pair_align',
            'Bioplus.permute_fasta',
            'Bioplus.prunetree',
            'Bioplus.sampler',
            'Bioplus.subset'
        ],
        version='0.1',
        author='Art Poon',
        url='https://github.com/PoonLab/Bioplus',
        classifiers=[
            'Programming Language :: Python :: 3',
            'License :: OSI Approved :: MIT License',
            'Operating System :: OS Independent'
        ]
    )


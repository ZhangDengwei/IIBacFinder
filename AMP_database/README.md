- overall.fa: combination of four databases
- overall.nr.fa: nonredundant sequences at 100% identity

# build index for DIAMOND
diamond makedb --in overall_nr.fa -d amp.index

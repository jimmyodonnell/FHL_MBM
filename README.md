# Methods in Marine Biodiversity

This is a collection of files for the Marine Biodiversity Methods course at Friday Habor Labs.

## Methods

### Taxonomic Annotation
We applied a number of quality filters to the sequence data (everything up to pre-clustering); now the aim is to apply a taxon name to each sequence.

Our approach for taxonomic annotation is the following:

- Does the sequence have a match in the FHL course barcode database at >97% identity? (blastn)
  - YES: use the name associated with that specimen
  - NO:
    - Does the sequence match a sequence in GenBank at >97% identity? (blastn)
      - YES: use the name associated with that GenBank sequence
      - NO:
        - After clustering at 5% (CROP), does OTU sequence have a match in GenBank at >85% identity? (blastn)
          - YES: use the phylum name associated with that GenBank sequence
          - NO:
            - Does the sequence have a probability of assignment to a sequence in (database?) with posterior probability >0.8?
              - YES: use the phylum name associated with that taxon
              - NO: sequence is labeled as 'unidentified taxon'

The output of this process is a table with the following columns; one row for each sequence in the pre-clustering fasta file:
  - seq_id [char]: sequence header from pre-clustering fasta file
  - FHL_db [boolean]: hit in local database?
  - GB_95 [boolean]: hit in GenBank >95%id?
  - CROP_ID [char]: sequence header of the OTU fasta file
  - GB_85 [boolean]: hit in GenBank >85%id after clustering?
  - SAP_prob [numeric,{0,1}]: posterior probability of assignment to a taxon using SAP
  - taxon_name [char]: taxon name 

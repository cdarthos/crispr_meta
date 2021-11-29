CREATE TABLE IF NOT EXISTS ncbi (
                Run_accession text NOT NULL UNIQUE,
                BioSample text,
                BioProject text,
                collection_date date,
                geo_loc_name text,
                isolate text,
                lat_lon text,
                sraID text,
                strain text,
                title text
                );
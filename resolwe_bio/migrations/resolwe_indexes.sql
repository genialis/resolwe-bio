-- Index for Django filter output__build.
CREATE INDEX idx_data_build ON flow_data ((output -> 'build'::text));

-- Index for Django filter output__feature_type.
CREATE INDEX idx_data_feature_type ON flow_data ((output -> 'feature_type'::text));

-- Index for Django filter output__source.
CREATE INDEX idx_data_source ON flow_data ((output -> 'source'::text));

-- Index for Django filter output__species__icontains.
CREATE INDEX idx_data_species_trgm ON flow_data USING gin (UPPER(output ->> 'species') gin_trgm_ops );

--
CREATE OR REPLACE FUNCTION generate_resolwe_bio_data_search(data_line flow_data)
    RETURNS tsvector
    LANGUAGE plpgsql
    AS $$
    DECLARE
        search tsvector;
    BEGIN
        SELECT
            -- Build.
            setweight(to_tsvector('simple', (output->>'build')), 'A') ||
            setweight(edge_ngrams((output->>'build')), 'B') ||
            setweight(edge_ngrams(get_characters((output->>'build'))), 'B') ||
            setweight(edge_ngrams(get_numbers((output->>'build'))), 'B') ||
            -- Feature type.
            setweight(to_tsvector('simple', (output->>'feature_type')), 'A') ||
            setweight(edge_ngrams((output->>'feature_type')), 'B') ||
            -- Source.
            setweight(to_tsvector('simple', (output->>'source')), 'A') ||
            setweight(edge_ngrams((output->>'source')), 'B') ||
            setweight(edge_ngrams(get_characters((output->>'source'))), 'B') ||
            setweight(edge_ngrams(get_numbers((output->>'source'))), 'B') ||
            -- Species.
            setweight(to_tsvector('simple', (output->>'species')), 'A') ||
            setweight(edge_ngrams((output->>'species')), 'B')
        FROM flow_data
        WHERE id=data_line.id
        INTO search;

        RETURN search;
    END;
    $$;

CREATE OR REPLACE FUNCTION data_biut()
    RETURNS TRIGGER
    LANGUAGE plpgsql
    AS $$
    BEGIN
        SELECT generate_resolwe_data_search(NEW) || generate_resolwe_bio_data_search(NEW)
        INTO NEW.search;

        RETURN NEW;
    END;
    $$;

-- Index for Django filter output__build.
CREATE INDEX idx_data_build ON flow_data ((output -> 'build'::text));

-- Index for Django filter output__feature_type.
CREATE INDEX idx_data_feature_type ON flow_data ((output -> 'feature_type'::text));

-- Index for Django filter output__source.
CREATE INDEX idx_data_source ON flow_data ((output -> 'source'::text));

-- Index for Django filter output__species__icontains.
CREATE INDEX idx_data_species_trgm ON flow_data USING gin (UPPER(output ->> 'species') gin_trgm_ops );


-- Add resolwe-bio specific fields to Data search.
CREATE OR REPLACE FUNCTION generate_resolwe_bio_data_search(data flow_data)
    RETURNS tsvector
    LANGUAGE plpgsql
    AS $$
    DECLARE
        search tsvector;
    BEGIN
        SELECT
            -- Build.
            setweight(to_tsvector('simple', COALESCE(data.output->>'build', '')), 'A') ||
            setweight(to_tsvector('simple', get_characters(data.output->>'build')), 'B') ||
            setweight(to_tsvector('simple', get_numbers(data.output->>'build')), 'B') ||
            -- Feature type.
            setweight(to_tsvector('simple', COALESCE(data.output->>'feature_type', '')), 'A') ||
            -- Source.
            setweight(to_tsvector('simple', COALESCE(data.output->>'source', '')), 'A') ||
            setweight(to_tsvector('simple', get_characters(data.output->>'source')), 'B') ||
            setweight(to_tsvector('simple', get_numbers(data.output->>'source')), 'B') ||
            -- Species.
            setweight(to_tsvector('simple', COALESCE(data.output->>'species', '')), 'A')
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

        IF NEW.search IS NULL THEN
            RAISE WARNING 'Search index for data (id: %) is NULL.', NEW.id;
        END IF;

        RETURN NEW;
    END;
    $$;

-- Add resolwe-bio specific fields to Entity search.
CREATE OR REPLACE FUNCTION generate_resolwe_bio_entity_search(entity flow_entity)
    RETURNS tsvector
    LANGUAGE plpgsql
    AS $$
    DECLARE
        search tsvector;
    BEGIN
        SELECT
            -- Species.
            setweight(to_tsvector('simple', COALESCE(entity.descriptor->'general'->>'species', '')), 'A')
        INTO search;

        RETURN search;
    END;
    $$;

CREATE OR REPLACE FUNCTION entity_biut()
    RETURNS TRIGGER
    LANGUAGE plpgsql
    AS $$
    BEGIN
        SELECT generate_resolwe_entity_search(NEW) || generate_resolwe_bio_entity_search(NEW)
        INTO NEW.search;

        IF NEW.search IS NULL THEN
            RAISE WARNING 'Search index for entity (id: %) is NULL.', NEW.id;
        END IF;

        RETURN NEW;
    END;
    $$;

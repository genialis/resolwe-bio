-- Trigger after insert/update Feature object.
CREATE OR REPLACE FUNCTION generate_resolwe_bio_kb_feature_search(feature resolwe_bio_kb_feature)
    RETURNS tsvector
    LANGUAGE plpgsql
    AS $$
    DECLARE
        search tsvector;
    BEGIN
        SELECT
            -- Feature name.
            setweight(to_tsvector('simple', feature.name), 'A') ||
            setweight(edge_ngrams(feature.name), 'B') ||
            -- Feature full name.
            setweight(to_tsvector('simple', feature.full_name), 'A') ||
            setweight(edge_ngrams(feature.full_name), 'B') ||
            -- Feature id.
            setweight(to_tsvector('simple', feature.feature_id), 'A') ||
            setweight(edge_ngrams(feature.feature_id), 'B') ||
            -- Feature aliases.
            setweight(to_tsvector('simple', array_to_string(feature.aliases, ' ')), 'B') ||
            setweight(edge_ngrams(array_to_string(feature.aliases, ' ')), 'C')
        INTO search;

        RETURN search;
    END;
    $$;

CREATE OR REPLACE FUNCTION resolwe_bio_kb_feature_biut()
    RETURNS TRIGGER
    LANGUAGE plpgsql
    AS $$
    BEGIN
        SELECT generate_resolwe_bio_kb_feature_search(NEW) INTO NEW.search;

        IF NEW.search IS NULL THEN
            RAISE WARNING 'Search index for kb feature (id: %) is NULL.', NEW.id;
        END IF;

        RETURN NEW;
    END;
    $$;

DROP TRIGGER IF EXISTS resolwe_bio_kb_feature_biut ON resolwe_bio_kb_feature;

CREATE TRIGGER resolwe_bio_kb_feature_biut
    BEFORE INSERT OR UPDATE
    ON resolwe_bio_kb_feature
    FOR EACH ROW EXECUTE PROCEDURE resolwe_bio_kb_feature_biut();

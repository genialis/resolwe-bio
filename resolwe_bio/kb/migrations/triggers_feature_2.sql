CREATE OR REPLACE FUNCTION generate_resolwe_bio_kb_feature_search(feature resolwe_bio_kb_feature)
    RETURNS tsvector
    LANGUAGE plpgsql
    AS $$
    DECLARE
        search tsvector;
    BEGIN
        SELECT
            -- Feature name.
            setweight(to_tsvector('simple_unaccent', feature.name), 'A') ||
            setweight(edge_ngrams(feature.name), 'B') ||
            -- Feature full name.
            setweight(to_tsvector('simple_unaccent', feature.full_name), 'B') ||
            setweight(edge_ngrams(feature.full_name), 'C') ||
            -- Feature id.
            setweight(to_tsvector('simple_unaccent', feature.feature_id), 'A') ||
            setweight(edge_ngrams(feature.feature_id), 'B') ||
            -- Feature aliases.
            setweight(to_tsvector('simple_unaccent', array_to_string(feature.aliases, ' ')), 'B') ||
            setweight(edge_ngrams(array_to_string(feature.aliases, ' ')), 'C')
        INTO search;

        RETURN search;
    END;
    $$;

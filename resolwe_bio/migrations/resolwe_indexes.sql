-- Index for Django filter output__build.
CREATE INDEX idx_data_build ON flow_data ((output -> 'build'::text));

-- Index for Django filter output__feature_type.
CREATE INDEX idx_data_feature_type ON flow_data ((output -> 'feature_type'::text));

-- Index for Django filter output__source.
CREATE INDEX idx_data_source ON flow_data ((output -> 'source'::text));

-- Index for Django filter output__species__icontains.
CREATE INDEX idx_data_species_trgm ON flow_data USING gin (UPPER(output ->> 'species') gin_trgm_ops );

-- Index for Django filter descriptor__subject_information__sample_label__icontains.
CREATE INDEX idx_entity_descriptor_sample_label_trgm ON flow_entity USING gin (UPPER(descriptor #>> '{subject_information,sample_label}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__subject_information__subject_id__icontains.
CREATE INDEX idx_entity_descriptor_subject_id_trgm ON flow_entity USING gin (UPPER(descriptor #>> '{subject_information,subject_id}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__subject_information__batch.
CREATE INDEX idx_entity_descriptor_batch ON flow_entity ((descriptor #>> '{subject_information,batch}'::text[]));

-- Index for Django filter descriptor__subject_information__group__iexact.
CREATE INDEX idx_entity_descriptor_group ON flow_entity (UPPER(descriptor #>> '{subject_information,group}'::text[]));

-- Index for Django filter descriptor__disease_information__disease_type__icontains.
CREATE INDEX idx_entity_descriptor_disease_type_trgm ON flow_entity USING gin (UPPER(descriptor #>> '{disease_information,disease_type}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__disease_information__disease_status__iexact.
CREATE INDEX idx_entity_descriptor_disease_status ON flow_entity (UPPER(descriptor #>> '{disease_information,disease_status}'::text[]));

-- Index for Django filter descriptor__immuno_oncology_treatment_type__io_drug__iexact.
CREATE INDEX idx_entity_descriptor_io_drug ON flow_entity (UPPER(descriptor #>> '{immuno_oncology_treatment_type,io_drug}'::text[]));

-- Index for Django filter descriptor__immuno_oncology_treatment_type__io_treatment__iexact.
CREATE INDEX idx_entity_descriptor_io_treatment ON flow_entity (UPPER(descriptor #>> '{immuno_oncology_treatment_type,io_treatment}'::text[]));

-- Index for Django filter descriptor__response_and_survival_analysis__confirmed_bor__iexact.
CREATE INDEX idx_entity_descriptor_confirmed_bor ON flow_entity (UPPER(descriptor #>> '{response_and_survival_analysis,confirmed_bor}'::text[]));

-- Index for Django filter descriptor__response_and_survival_analysis__pfs_event__iexact.
CREATE INDEX idx_entity_descriptor_pfs_event ON flow_entity (UPPER(descriptor #>> '{response_and_survival_analysis,pfs_event}'::text[]));

-- Index for Django filter descriptor__general__description__icontains.
CREATE INDEX idx_entity_descriptor_description_trgm ON flow_entity USING gin (UPPER(descriptor #>> '{general,description}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__general__biosample_source__icontains.
CREATE INDEX idx_entity_descriptor_biosample_source_trgm ON flow_entity USING gin (UPPER(descriptor #>> '{general,biosample_source}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__general__biosample_treatment__icontains.
CREATE INDEX idx_entity_descriptor_biosample_treatment_trgm ON flow_entity USING gin (UPPER(descriptor #>> '{general,biosample_treatment}'::text[]) gin_trgm_ops );

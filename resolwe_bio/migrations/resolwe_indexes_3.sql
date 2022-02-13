-- Index for Django filter descriptor__general__species__icontains
CREATE INDEX idx_entity_descriptor__general__species ON flow_entity USING gin (UPPER(descriptor #>> '{general,species}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__general__organ__icontains
CREATE INDEX idx_entity_descriptor__general__organ ON flow_entity USING gin (UPPER(descriptor #>> '{general,organ}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__general__biosample_source__icontains
CREATE INDEX idx_entity_descriptor__general__biosample_source ON flow_entity USING gin (UPPER(descriptor #>> '{general,biosample_source}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__disease_information__organ_part__icontains
CREATE INDEX idx_entity_descriptor__disease_information__organ_part ON flow_entity USING gin (UPPER(descriptor #>> '{disease_information,organ_part}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__disease_information__biopsy_site__icontains
CREATE INDEX idx_entity_descriptor__disease_information__biopsy_site ON flow_entity USING gin (UPPER(descriptor #>> '{disease_information,biopsy_site}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__pathological_information__organ_part__icontains
CREATE INDEX idx_entity_descriptor__pathological_information__organ_part ON flow_entity USING gin (UPPER(descriptor #>> '{pathological_information,organ_part}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__pathological_information__biopsy_site__icontains
CREATE INDEX idx_entity_descriptor__pathological_information__biopsy_site ON flow_entity USING gin (UPPER(descriptor #>> '{pathological_information,biopsy_site}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__general__biosample_treatment__icontains
CREATE INDEX idx_entity_descriptor__general__biosample_treatment ON flow_entity USING gin (UPPER(descriptor #>> '{general,biosample_treatment}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__treatment_type__drug__icontains
CREATE INDEX idx_entity_descriptor__treatment_type__drug ON flow_entity USING gin (UPPER(descriptor #>> '{treatment_type,drug}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__immuno_oncology_treatment_type__io_drug__icontains
CREATE INDEX idx_entity_descriptor__immuno_oncology_treatment_type__io_drug ON flow_entity USING gin (UPPER(descriptor #>> '{immuno_oncology_treatment_type,io_drug}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__response_and_survival_analysis__clinical_benefit__icontains
CREATE INDEX idx_entity_descriptor__response_and_survival_analysis__clinical_benefit ON flow_entity USING gin (UPPER(descriptor #>> '{response_and_survival_analysis,clinical_benefit}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__response_and_survival_analysis__confirmed_bor__icontains
CREATE INDEX idx_entity_descriptor__response_and_survival_analysis__confirmed_bor ON flow_entity USING gin (UPPER(descriptor #>> '{response_and_survival_analysis,confirmed_bor}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__response_and_survival_analysis__unconfirmed_bor__icontains
CREATE INDEX idx_entity_descriptor__response_and_survival_analysis__unconfirmed_bor ON flow_entity USING gin (UPPER(descriptor #>> '{response_and_survival_analysis,unconfirmed_bor}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__response_and_survival_analysis__pfs__icontains
CREATE INDEX idx_entity_descriptor__response_and_survival_analysis__pfs ON flow_entity USING gin (UPPER(descriptor #>> '{response_and_survival_analysis,pfs}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__response_and_survival_analysis__os__icontains
CREATE INDEX idx_entity_descriptor__response_and_survival_analysis__os ON flow_entity USING gin (UPPER(descriptor #>> '{response_and_survival_analysis,os}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__response_and_survival_analysis__dfs__icontains
CREATE INDEX idx_entity_descriptor__response_and_survival_analysis__dfs ON flow_entity USING gin (UPPER(descriptor #>> '{response_and_survival_analysis,dfs}'::text[]) gin_trgm_ops );

-- Index for Django filter descriptor__response_and_survival_analysis__ttp__icontains
CREATE INDEX idx_entity_descriptor__response_and_survival_analysis__ttp ON flow_entity USING gin (UPPER(descriptor #>> '{response_and_survival_analysis,ttp}'::text[]) gin_trgm_ops );

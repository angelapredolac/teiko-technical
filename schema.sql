PRAGMA foreign_keys = ON;

CREATE TABLE projects (
    project TEXT PRIMARY KEY NOT NULL
);
CREATE TABLE subjects (
    subject TEXT PRIMARY KEY NOT NULL,
    project TEXT NOT NULL REFERENCES projects(project),
    condition TEXT NOT NULL,
    age INTEGER NOT NULL CHECK(age >= 0),
    sex TEXT NOT NULL CHECK(sex IN ('M', 'F')),
    treatment TEXT NOT NULL,
    response TEXT CHECK(response IN ('yes', 'no'))
);
CREATE TABLE samples (
    sample TEXT PRIMARY KEY NOT NULL,
    subject TEXT NOT NULL REFERENCES subjects(subject),
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL
);
CREATE TABLE populations (
    population TEXT PRIMARY KEY NOT NULL
);
CREATE TABLE cell_counts (
    sample TEXT NOT NULL REFERENCES samples(sample),
    population TEXT NOT NULL REFERENCES populations(population),
    count INTEGER NOT NULL CHECK(count >= 0),
    PRIMARY KEY(sample, population)
);
CREATE INDEX samples_subject_idx ON samples(subject);
CREATE INDEX subjects_cohort_idx ON subjects(condition, treatment);

CREATE VIEW sample_frequencies AS
WITH totals AS (
    SELECT sample, SUM(count) AS total_count FROM cell_counts GROUP BY sample
)
SELECT c.sample, t.total_count, c.population, c.count,
       100.0 * c.count / NULLIF(t.total_count, 0) AS percentage
FROM cell_counts c JOIN totals t USING(sample);

CREATE VIEW sample_metadata AS
SELECT s.*, u.project, u.condition, u.age, u.sex, u.treatment, u.response
FROM samples s JOIN subjects u USING(subject);

CREATE VIEW baseline_samples AS
SELECT * FROM sample_metadata
WHERE condition = 'melanoma' AND treatment = 'miraclib'
  AND sample_type = 'PBMC' AND time_from_treatment_start = 0;

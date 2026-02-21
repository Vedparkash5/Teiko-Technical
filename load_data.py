import sqlite3
import csv
import os

#set the root directory path
CSV_PATH = "cell-count.csv"
DB_PATH = "cell_counts.db"

# connect the database
conn = sqlite3.connect(DB_PATH)
cursor = conn.cursor()

# -----------  Database Schema  -----------
# 2 tables, subjects and samples

# Entities
## Subject : subject (primary key), project, condition, age, sex, treatment, response
## Sample : sample (primary Key), subject(foreign key), sample_type, time_from_treatment_start, 5 cells

#create the database
cursor.executescript("""
    CREATE TABLE IF NOT EXISTS subjects (
        subject     TEXT PRIMARY KEY,
        project     TEXT,
        condition   TEXT,   
        age         INTEGER,    
        sex         TEXT,   
        treatment   TEXT,
        response    TEXT
    );

    CREATE TABLE IF NOT EXISTS samples (
        sample                      TEXT PRIMARY KEY,
        subject                     TEXT,
        sample_type                 TEXT,
        time_from_treatment_start   INTEGER,
        b_cell                      INTEGER,
        cd8_t_cell                  INTEGER,
        cd4_t_cell                  INTEGER,
        nk_cell                     INTEGER,
        monocyte                    INTEGER,
        FOREIGN KEY (subject) REFERENCES subjects(subject)
    );
""")

with open(CSV_PATH, newline="") as f:
    reader = csv.DictReader(f)

    subjects_seen = set()

    for row in reader:
        subject_id = row["subject"]

        #insert into subjects (only once per subject)
        if subject_id not in subjects_seen:
            cursor.execute("""
                INSERT OR IGNORE INTO subjects
                (subject, project, condition, age, sex, treatment,response)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (
                subject_id,
                row["project"],
                row["condition"],
                int(row["age"]),
                row["sex"],
                row["treatment"],
                row["response"]
            ))
            subjects_seen.add(subject_id)

        cursor.execute("""
                INSERT OR IGNORE INTO samples
                (sample, subject, sample_type, time_from_treatment_start, b_cell, cd8_t_cell, cd4_t_cell, nk_cell,monocyte)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                row["sample"],
                subject_id,
                row["sample_type"],
                int(row["time_from_treatment_start"]),
                int(row["b_cell"]),
                int(row["cd8_t_cell"]),
                int(row["cd4_t_cell"]),
                int(row["nk_cell"]),
                int(row["monocyte"])
            ))

conn.commit()
conn.close()
print("Database loaded successfully!")
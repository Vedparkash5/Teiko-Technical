# Teiko-Technical

# Clinical Trial — Immune Cell Population Analysis For Medicine Effectiveness

By Vedparkash Singh  

---

## Quickstart

### GitHub Codespaces 

1. Open Github Codespaces from the repo
2. The necassary files should auto download.
3. The database will be created by running `python3 load_data.py`
4. Run the dashboard:
   ```bash
   python -m streamlit run dashboard.py
   ```
   You should see a popup on the bottom right corner which says 'Your application running on port 8501 is available'. Click the green button that says open in browser. Dashboard will open on a new page.
5. Analysis code for Parts 2-4 can be found on analysis.ipynb Jupyter notebook 

---

## Database Schema

I decided to create relational database with 2 tables. The csv file has all of the data pertaining to the study, which includes experimental methodologies, subject conditions, subject responses, subject cell counts, and subject demographics. I decided to create two tables, one for data involving the subjects and another for data pertaining to the sampled experiment. 

In terms of scalability, this system allows the scientists (Bob and Yah) to perform a lot more trials. As trail number grows, new sample data can be added to the samples table without duplicating records into the subjects table. Furthermore, if we want to expand testing to different parameters and measurements, we can add them to our existing tables by creating new columns, both for subjects and samples. For example, if we want to add additional traits such as activity level, weight, etc to each subjects records, they can easily just be appended as new columns. 


---

## Analysis Overview

### Part 1 — Data Loading (`load_data.py`)
Reads `cell-count.csv` and loads it into a normalized SQLite database. Each subject is inserted once into `subjects`; each sample row is inserted into `samples` with a foreign key reference.

### Part 2 — Frequency Table (`analysis.py`)
First, I extracted the data from the Database and created two dataframes for each table. I added a total_count column for each cell type for each sample. Then for each sample, I calculated the relative frequency of each immune cell population as a percentage of the total cell count across all five populations. Output is a table with one row per population per sample.

### Part 3 — Statistical Analysis 
I compared the cell population relative frequencies between responders and non-responders in the melanoma + miraclib + PBMC cohort.<br>

**Dataframe Creation**<br>
I joined the 2 tables by subject and filtered it by the constraints set by Bob. Then created a new datafame with the filtered data.<br>

I compared the mean and variance statics between responders and non-responders. To further this exploration, I created a boxplot comparing the data distribution between responders and non-responders for each cell type.<br>

**Statistical Tests:**<br>
First, in order to determine what statistical tests to perform, I analyzed the distribution of data for all data pools using the Shapiro-Wilk Normality Test. The data is mostly not normally distributed. <br>

I chose to perform 3 analytics tests that work great for data that is not uniformally distributed between 2 independent groups. Mann-Whitney and Brunner-Munzel tests are perfect for this. I also did a Welch's T test t because it evaluates the actual values, not just distributions.<br>

Additionally, because there are 5 cell counts that we are testing, we have to apply the Bonferroni correction to adjust the threshold from 0.05 to 0.01 for 5 simultaneous comparisions.<br>

The key conclusions I found is that the results are unconclusive with the current data pool, non of the tests passed under the Bonferroni threshold. The cd4_t cell was the only cell population that came close to the threshold with values 0f 0.013 for the Mann-Whitney and Brunner-Munzel tests, it has a passing score for the Welch's T test. This warrants further investigation. cd4_t cell also had the largest average population so maybe more indept research needs to be performed on different aspects of the cell type in patients, especially relating to the effectiveness of the medicines.<br>


### Part 4 — Baseline Subset Tests<br>
Analysis restricted to baseline samples was very straight forward. I extracted the data from the database by joining the 2 tables and filtering out the specifications of melanoma patients that were treated with miraclib at time = 0. Once the dataframe was made, for each request for data on samples per project, subjects per response, and subjects per gender was a simple grouping and counting function.<br>

### Dashboard (`dashboard.py`)<br>
I built the dashboard as a single-page scrollable Streamlit report organized as a knowledge journey through seven sections: The Data → Cell Population Overview → Responders vs. Non-Responders → Statistical Analysis → Baseline Cohort → Key Findings → Future Directions. I built the dashboard by mixing computation with markup html code chunks in order to make the dashboard informative yet a little asthetic. I still kept the styling very professional similar to a Jupyter notebook. 




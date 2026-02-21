# Teiko-Technical

# Clinical Trial — Immune Cell Population Analysis For Medicine Effectiveness

By Vedparkash Singh  

---

## Quickstart

### GitHub Codespaces 

1. Open Github Codespaces from the repo 
4. The database will be created by running `python3 load_data.py`
5. Run the dashboard:
   ```bash
   python -m streamlit run dashboard.py
   ```
6. Analysis code for Parts 2-4 can be found on analysis.ipynb Jupyter notebook 

---

## Database Schema

I decided to create relational database with 2 tables. The csv file has all of the data pertaining to the study, which includes experimental methodologies, subject conditions, subject responses, subject cell counts, and subject demographics. I decided to create two tables, one for data involving the subjects and another for data pertaining to the sampled experiment. 

In terms of scalability, this system allows the scientists (Bob and Yah) to perform a lot more trials. As trail number grows, new sample data can be added to the samples table without duplicating records into the subjects table. Furthermore, if we want to expand testing to different parameters and measurements, we can add them to our existing tables by creating new columns, both for subjects and samples. For example, if we want to add additional traits such as activity level, weight, etc to each subjects records, they can easily just be appended as new columns. 



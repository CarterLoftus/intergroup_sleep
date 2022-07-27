# intergroup_sleep

This repository contains all the code needed to reproduce the results presented in Chapter 3 of my dissertation (The Physiology, Behavioral Ecology, and Collective Dynamics of Sleep in Wild Olive Baboons (Papio anubis)). The README presented here explains where the raw data is located and provides the instructions needed to repeat the data processing and analysis, allowing a reader to completely reproduce the results.

All analyses are based on the GPS and accelerometry data contained within the Movebank project named “Leopards, vervets, and baboons in Laikipia, Kenya”. Please contact Lynne Isbell (laisbell@ucdavis.edu) for permission and access to download the data. The analyses also involved observations of animations of the GPS data during inter-group encounters (see below). These observations are available for download at Dryad (https://doi.org/10.25338/B87D15).

To repeat the analyses, clone this repository to a local computer, and copy the data into a directory titled "DATA" within the repository. Then run the scripts in the order they are numbered.


The script “00_data_cleaning.R” loads the GPS data from Movebank, converts the longitude-latitude coordinates to UTM coordinates, separates the data by species, and cleans the data. The data cleaning consists of interpolating GPS fixes to the ideal (i.e. programmed) sampling schedule, interpolating missing location data, and replacing anomalous location data with interpolated locations. The cleaned data is then output as “bab_df.csv”

The script “01_home_ranges_and_overlap.R” uses “bab_df.csv” as an input and calculates the groups’ home ranges as well as the overlap and Bhattacharyya's affinity between them with the “adehabitatHR” package in R. 

The script “02_individual_movement_properties_and_sleep_sites.R” takes “bab_df.csv” as an input and calculates some basic properties of movement and identifies sleep sites from the GPS data. Specifically, it discretizes the movement tracks of each individual to a 30 m resolution, and in doing so, calculates spatially discretized step lengths. This script also calculates headings and turning angles of each individual. It also identifies distinct sleep sites, and determines the distinct sleep site in which each individual spent each night. All of this information is then output in the data frame, “bab_complete.csv”.  The script also identifies instances of group co-sleeping (i.e. when two groups share a sleep site).

The script “03_prepping_ACC_data_for_sleep_analysis.R” takes the full dataset downloaded directly from Movebank as an input (download data with all sensor data from the project “Leopards, vervets, and baboons in Laikipia, Kenya”, and stored in the “DATA” folder within repository. The input csv should be named, by Movebank’s default, “Leopards, vervets, and baboons in Laikipia, Kenya.csv”). Using this input, the script calculates the average VeDBA (vectorial dynamic body acceleration) over each accelerometry burst. The script then takes the log of the average VeDBA for each burst, and stores this information in the output data frame “full_night_and_day_data.csv”. 

The script “04_effect_of_cosleeping_on_sleep_quality.R” extracts sleep metrics and tests the effects of sleeping with another group on total sleep time, sleep efficiency, and sleep fragmentation. This script also plots the model predictions that are represented in Figure 3.1 of the dissertation.

The script “05_effect_of_cosleeping_on_sleep_synchronization.R” takes the minute-by-minute sleep classification (“full_dat.csv”), as well as “bab_complete.csv” as an input, and models the influence of two groups sleeping in the same sleep site on the amount of synchronization in the sleep-wake patterns during the night between non-group-mates. The script also plots the predictions of the model, although the plot does not appear as a figure in the dissertation.

The script “06_dyadic_movement_properties.R” takes “bab_complete.csv” as an input and calculates the distance between each dyad of individuals at each timestamp, as well as the dynamic time warping distance between the trajectories of individuals with a sliding window approach. It stores this information in the output data frame “bab_dyad.csv”.

The script “07_identifying_response_distance.R” takes “bab_complete.csv” as an input and determines the distance that baboon groups respond to neighboring groups, using a permutation technique. It stores the results of these permutations in “emp_detect_dist_100m_7000_seed_500.csv” and “perm_detect_dist_100m_7000_seed_500.csv”, and these data frames are used within this script to identify the response distance (which was determined to be 600m).

The script “08_extracting_encounters.R” takes “bab_dyad.csv” as an input and uses it to extract encounters that happen between study groups, based on the 600m encounter threshold determined in “07_identifying_response_distance.R”. This script creates the data frame “baboon_dyadDurat.csv”, which gives a complete account of all the inter-group encounters (with the identities of the groups involved, as well as the start time, end time, and duration of encounters). The script then goes on to use this “baboon_dyadDurat.csv” data frame, along with the GPS data from “baboon_complete.csv” to produce KMLs of locations of all study individuals during each encounter.

I used the KMLs produced in “08_extracting_encounters.R” to then score the interactions that occurred within each encounter. I manually added columns that contained the manual scorings for these interactions to the “baboon_dyadDurat.csv”, and saved the manually extended data frame as “dyadDurat_interaction_scoring.csv”. This data frame is available on Dryad (https://doi.org/10.25338/B87D15).

The script “09_parsing_interaction_scoring.R”  takes “bab_complete.csv” and “dyadDurat_interaction_scoring.csv” (see previous paragraph) as inputs, and parses the interaction scoring into a format that allowed us to model the data. This parsed data is then saved as “dyad_interact_df.csv”.

The script “10_effect_of_cosleeping_on_next_day_interactions.R” takes “bab_complete.csv” and “dyad_interact_df.csv” as inputs and prepares a data frame for modelling. The script then models the influence of two groups sharing the same sleep site on the probability that they will interact the following day and on the probability that, if they interact, the interaction will involve cohesive movement. The script also plots the model predictions that are shown in Figure 3.2 in the dissertation.

The script “11_movement_path_randomizations.R” takes “bab_complete.csv” as an input and performs the movement path permutation procedure that we used to test whether baboon groups’ movements with respect to each other deviated from that expected by random chance. The script both performs the permutations as well as evaluates significance and plots results (Figure 3.3 in dissertation).

The script “12_movement_path_randomizations_cosleep_split.R” takes “bab_complete.csv” as an input and performs the movement path permutation procedures as in “11_movement_path_randomizations.R”, but this script analyses separately the days following nights of group co-sleeping and days following nights on which groups slept as lone groups. The script performs the permutations, evaluates significance (for both days following co-sleeping and days following sleeping alone), and plots the results.

The script “13_sleep_site_occupancy_randomizations” takes “bab_complete.csv” as an input, and performs the permutations that we used to test whether and how baboon groups’ frequency of sleep site sharing deviated from that expected by chance. The script also evaluates significance and plots the results. 





# FastSimNN
Instant simulation of detector responses using ML.
Unlike the traditional Geant4 and fast detector simulation programs,
this approach makes almost instant transformation (1-2 ms/event) of true-level Monte Carlo
record to a record that closely follows the full/fast detector simulation.
It uses an array of neural networks that need to be trained beforehand using true-level and
detector-level records.


# Now to run
The program is designed and compiled to run on the ANL ATLAS cluster.
The input of this program are data from the HepSim repository.


Version 1.0, December 1, 2018: Initial version
Version 2.0, December 14, 2018: Corrected Eta resolution

Send  comments to S.Chekanov (ANL)

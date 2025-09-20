# The Integrated Healthcare Timetabling Competition

This project focuses on the challenge of integrated healthcare scheduling, an issue in hospital management where multiple scheduling decisions must be coordinated to ensure efficient patient care. It is based on the
problem description and rules from the [Integrated Healthcare Timetabling Competition 2024](https://ihtc2024.github.io/), with adaptations to suit the specific needs of this project.

The project combines two scheduling challenges:
1. Patient Admission Scheduling – Determining when patients are admitted and assigning them to hospital rooms.
2. Nurse-to-Room Assignment – Scheduling nurses to cover patient care in different hospital rooms.

# Solution

A Mixed-Integer Linear Program (MILP) was formulated mathematically. Then, it was solved using the Gurobi solver, using Julia as an interface. Heuristics were also implemented.

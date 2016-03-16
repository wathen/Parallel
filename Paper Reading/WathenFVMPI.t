1. What problem does the paper address?

The paper addresses the problem of formal verification of MPI programs. They do this without building verification models.

2. What is the key insight/idea in the paper's solution to this problem?

The key insight is to use ISP to create a "push-button" dynamic model checker.

3. What did the authors do to demonstrate their claims? (e.g. implement a tool, present a proof, run simulations, etc.)

The authors demonstrate there claims by creating a "push-button" dynamic model checker for MPI programs. The authors claim that this is the first formal verification tool that handles C/MPI programs.

4. Is the support for the claims convincing? Why or why not?

They create a formal verification too for two case studies: ParMETIS and MADRE.There results show that these case studies are free from deadlocks and other breakdowns.

5. What are your questions or other comments about the paper?

Having tools that can verify MPI codes is very useful, especially for a scientific computing person. The majority of our parallel codes are based on an MPI framework. For large problems, in excess of 100 million degrees of freedom, we cannot afford our MPI based code to stall and crash. Verification is therefore very important.

```mermaid
graph TD;
  A(Start) --> B{"Does your data contain 0 (CT = 40) values?"};
  B -->|Yes| C{Do your data represent single-digit copy numbers?};
  B -->|No| D(`Use a standard model from another package (e.g. glm(), lm(), etc.)`);
  C -->|Yes| E(`Use the Poisson-distributed GLM models (eDNA_count)`);
  C -->|No| F{Do the zeros in the data occur in otherwise high concentrations?"};
  F -->|Yes| G(Use a zero-inflated model (eDNA_lm*_zinf()));
  F -->|No| H(Use a censored-data model (eDNA_lm*()));
  
```
	

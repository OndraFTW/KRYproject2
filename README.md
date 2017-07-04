# KRYproject2
Generování klíčů, dešifrování, šifrování a prolomení klíčů algoritmu RSA.

## Použití
### Generování klíčů
vstup: ./kry -g B  
výstup: P Q N E D

### Šifrování
vstup: ./kry -e E N M  
výstup: C

### Dešifrování
vstup: ./kry -d D N C  
výstup: M

### Prolomení RSA
vstup: ./kry -b N  
výstup: P

### Legenda
B ... požadovaná velikost veřejného modulu v bitech (např. 1024). Nebude jej tedy možné zapsat na méně jak B bitů.  
P ... prvočíslo (při generování náhodné)  
Q ... jiné prvočíslo (při generování náhodné)  
N ... veřejný modulus  
E ... veřejný exponent (většinou 3)  
D ... soukromý exponent  
M ... otevřená zpráva (číslo, nikoli text)  
C ... zašifrovaná zpráva (číslo, nikoli text)  

Všechna čísla na vstupu i výstupu (kromě B) jsou hexadecimální, začínají "0x" výstup končí novým řádkem.

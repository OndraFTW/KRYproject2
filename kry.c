/*
	Autor: Ondřej Šlampa, xslamp01@stud.fit.vutbr.cz
	Projekt: KRY 2
	Popis: Program pro práci s RSA šifrováním.
*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<stdbool.h>
#include<errno.h>
#include<time.h>
#include<gmp.h>

/*
	Akce zadaná na příkazové řádce.
*/
typedef enum action{
    ERROR,GENERATE,ENCRYPT,DECRYPT,BREAK
} Action;

/*
	Argumenty příkazové řádky.
*/
typedef struct args{
    Action a;
    char* b;
    char* key;
    char* n;
    char* msg;
} Args;

/*
	Zpracování argumentů příkazové řádky.
	argc počet argumentů
	argv argumenty
*/
Args get_args(int argc, char** argv){
    Args args;
    //generování klíčů
    if(argc==3 && strcmp(argv[1], "-g")==0){
        args.a=GENERATE;
        args.b=argv[2];
    }
    //šifrování
    else if(argc==5 && strcmp(argv[1], "-e")==0){
        args.a=ENCRYPT;
        args.key=argv[2];
        args.n=argv[3];
        args.msg=argv[4];
    }
    //dešifrování
    else if(argc==5 && strcmp(argv[1], "-d")==0){
        args.a=DECRYPT;
        args.key=argv[2];
        args.n=argv[3];
        args.msg=argv[4];
    }
    //prolomení
    else if(argc==3 && strcmp(argv[1], "-b")==0){
        args.a=BREAK;
        args.n=argv[2];
    }
    else{
        args.a=ERROR;
    }
    return args;
}

/*
	Vrátí počáteční hodnotu pro generátor náhodných čísel.
*/
unsigned long seed(){
    return clock();
}

/*
	Nalezne inverzní prvek ve zbytkové třídě poimocí rozšířeného Eukleidova
	algoritmu.
	ret inverzní prvek prvku a
	a prvek
	n modulo zbytkové třídy
*/
void inverse(mpz_t ret, mpz_t a, mpz_t n){
    mpz_t t;
    mpz_init(t);
    mpz_set_ui(t,0ul);
    
    mpz_t newt;
    mpz_init(newt);
    mpz_set_ui(newt,1ul);
    
    mpz_t r;
    mpz_init(r);
    mpz_set(r,n);
    
    mpz_t newr;
    mpz_init(newr);
    mpz_set(newr,a);
    
    mpz_t quotient;
    mpz_init(quotient);
    
    mpz_t tmp;
    mpz_init(tmp);
    
    while(mpz_cmp_ui(newr,0ul)!=0){
        mpz_fdiv_q(quotient,r,newr);
        
        mpz_set(tmp,t);
        mpz_set(t,newt);
        mpz_submul(tmp,quotient,newt);
        mpz_set(newt,tmp);
        
        mpz_set(tmp,r);
        mpz_set(r,newr);
        mpz_submul(tmp,quotient,newr);
        mpz_set(newr,tmp);
    }
    
    if(mpz_cmp_ui(r,1ul)>0){
        mpz_set_ui(ret,0ul);
    }
    else if(mpz_cmp_ui(t,0ul)<0){
        mpz_add(ret,t,n);
    }
    else{
        mpz_set(ret,t);
    }
    
    mpz_clear(t);
    mpz_clear(newt);
    mpz_clear(r);
    mpz_clear(newr);
    mpz_clear(quotient);
    mpz_clear(tmp);
}

/*
	Nalezne jacobiho symbol pro (a_in/n_in) ve zbytkové třídě mod.
	Funguje pro liché n && a<n.
	r výsledný jacobiho symbol
	a_in čitatel
	n_in jmenovatel
	mod modulo zbytkové třídy
*/
void jacobi(mpz_t r, mpz_t a_in, mpz_t n_in, mpz_t mod){
    int j=1;
    
    mpz_t a;
    mpz_init(a);
    mpz_set(a,a_in);
    
    mpz_t n;
    mpz_init(n);
    mpz_set(n,n_in);
    
    mpz_t tmp;
    mpz_init(tmp);
    
    while(mpz_cmp_ui(a, 0ul)!=0){
        while(mpz_even_p(a)){
            mpz_fdiv_q_ui(a,a,2ul);
            mpz_mod_ui(tmp,n,8ul);
            if(mpz_cmp_ui(tmp, 3ul)==0 || mpz_cmp_ui(tmp, 5ul)==0){
                j=-j;
            }
        }
        mpz_swap(a,n);
        mpz_mod_ui(tmp,a,4ul);
        if(mpz_cmp_ui(tmp, 3ul)==0){
            mpz_mod_ui(tmp,n,4ul);
            if(mpz_cmp_ui(tmp, 3ul)==0){
                j=-j;
            }
        }
        mpz_mod(a,a,n);
    }
    
    if(mpz_cmp_ui(n,1ul)==0){
        if(j==1){
            mpz_set_ui(r,1ul);
        }
        else{
            mpz_sub_ui(r,mod,1ul);
        }
    }
    else{
        mpz_set_ui(r,0ul);
    }
    
    mpz_clear(a);
    mpz_clear(n);
    mpz_clear(tmp);
}

/*
	Zjistí jesli je n prvočíslo. Při odpovědi true je šance 1:2100, že n není
	prvočíslo.
	n vstup
*/
bool is_prime(mpz_t n){
    
    if(mpz_even_p(n)){
        return false;
    }
    
    if(mpz_cmp_ui(n,0ul)<0){
        return false;
    }
    
    if(mpz_cmp_ui(n,3ul)==0){
        return true;
    }
    
    mpz_t a;
    mpz_init(a);
    
    mpz_t x;
    mpz_init(x);
    
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_t top;
    mpz_init(top);
    mpz_sub_ui(top,n,3ul);

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,seed());
    
    bool r=true;
    
    //gmp_printf("n=%Zd\n",n);
    for(int i=0;i<100;i++){
        mpz_urandomm(a,state,top);
        mpz_add_ui(a,a,2ul);
        //gmp_printf("a=%Zd\n",a);
        
        jacobi(x,a,n,n);
        //gmp_printf("x=%Zd\n",x);
        
        if(mpz_cmp_ui(x, 0ul)==0){
            r=false;
            break;
        }
        
        mpz_sub_ui(tmp,n,1ul);
        mpz_fdiv_q_ui(tmp,tmp,2ul);
        mpz_powm(tmp,a,tmp,n);
        
        if(mpz_cmp(tmp,x)!=0){
            r=false;
            break;
        }
    }
    
    mpz_clear(a);
    mpz_clear(x);
    mpz_clear(tmp);
    mpz_clear(top);
    gmp_randclear(state);
    
    return r;
}

/*
	Vygeneruje prvočíslo n l bitech.
	p výsledné prvočíslo
	l počet bitů prvočísla
*/
void get_prime(mpz_t p, unsigned long l){
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,seed());
    
    do{
        mpz_urandomb(p,state,(mp_bitcnt_t)l);
    }while(!is_prime(p));
    
    gmp_randclear(state);
}

/*
	Nalezení největšího společného dělitele.
	r výsledný největší společný dělitel u_in a v_in
	u_in číslo
	v_in číslo
*/
void gcd(mpz_t r, mpz_t u_in, mpz_t v_in){
	unsigned long shift;
	
	if(mpz_cmp(u_in,v_in)==0){
		mpz_set(r,u_in);
		return;
	}
	
	if(mpz_cmp_ui(u_in,0ul)==0){
		mpz_set(r,v_in);
		return;
	}
	
	if(mpz_cmp_ui(v_in,0ul)==0){
		mpz_set(r,u_in);
		return;
	}
	
	mpz_t u;
    mpz_init(u);
    mpz_set(u, u_in);
	
	mpz_t v;
    mpz_init(v);
    mpz_set(v, v_in);
	
	mpz_t tmp;
    mpz_init(tmp);
	
	mpz_t one;
    mpz_init(one);
    mpz_set_ui(one,1ul);
	
	for(shift=0;mpz_cmp_ui((mpz_ior(tmp,u,v),mpz_and(tmp,tmp,one),tmp),0ul)==0;++shift){
		mpz_tdiv_q_2exp(u,u,1ul);
		mpz_tdiv_q_2exp(v,v,1ul);
	}
	
	while(mpz_cmp_ui((mpz_and(tmp,u,one),tmp),0ul)==0){
		mpz_tdiv_q_2exp(u,u,1ul);
	}
	
	do{
		while(mpz_cmp_ui((mpz_and(tmp,v,one),tmp),0ul)==0){
			mpz_tdiv_q_2exp(v,v,1ul);
		}
		
		if(mpz_cmp(u,v)>0){
			mpz_swap(u,v);
		}
		
		mpz_sub(v,v,u);
	}
	while(mpz_cmp_ui(v,0ul)!=0);
	
	mpz_mul_2exp(r,u,shift);
	
	mpz_clear(u);
	mpz_clear(v);
	mpz_clear(tmp);
	mpz_clear(one);
}

/*
	Generování RSA klíčů.
	Vrací EXIT_SUCCESS nebo EXIT_FAILURE.
	args argumenty příkazové řádky
*/
int generate(Args args){
    char* endptr=NULL;
    unsigned long l=strtoul(args.b,&endptr,10);
    
    if(errno==ERANGE || strcmp(endptr,"")!=0){
        return EXIT_FAILURE;
    }
    
    mpz_t p;
    mpz_init(p);
    
    mpz_t q;
    mpz_init(q);
    
    mpz_t n;
    mpz_init2(n,l);
    
    do{
    	get_prime(p,l/2+1);
    	get_prime(q,l/2+1);
    	mpz_mul(n,p,q);
    }
    while(mpz_sizeinbase(n,2)!=l);
    
    mpz_t phi_n;
    mpz_init(phi_n);
    mpz_set(phi_n,n);
    mpz_sub(phi_n,phi_n,p);
    mpz_sub(phi_n,phi_n,q);
    mpz_add_ui(phi_n,phi_n,1ul);
    
    mpz_t e;
    mpz_init(e);
    
    mpz_t tmp;
    mpz_init(tmp);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,seed());
    
    do{
    	mpz_urandomm(e,state,phi_n);
    }
    while(mpz_cmp_ui((gcd(tmp,e,phi_n),tmp),1ul)!=0);
    
    mpz_t d;
    mpz_init(d);
    inverse(d,e,phi_n);
    
    //gmp_printf("phi=%#Zx\n",phi_n);
    gmp_printf("%#Zx %#Zx %#Zx %#Zx %#Zx\n",p,q,n,e,d);
    
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(n);
    mpz_clear(phi_n);
    mpz_clear(e);
    mpz_clear(d);
    mpz_clear(tmp);
    gmp_randclear(state);
    
    return EXIT_SUCCESS;
}

/*
	Šifrování a dešifrování pomocí RSA.
	Vrací EXIT_SUCCESS nebo EXIT_FAILURE.
	args argumenty příkazové řádky
*/
int encrypt(Args args){
    int ctrl=-1;
    
    mpz_t key;
    mpz_init(key);
    ctrl=mpz_set_str(key,args.key,0);
    if(ctrl!=0){
        mpz_clear(key);
        return EXIT_FAILURE;
    }
    
    mpz_t n;
    mpz_init(n);
    ctrl=mpz_set_str(n,args.n,0);
    if(ctrl!=0){
        mpz_clear(key);
        mpz_clear(n);
        return EXIT_FAILURE;
    }
    
    mpz_t msg;
    mpz_init(msg);
    ctrl=mpz_set_str(msg,args.msg,0);
    if(ctrl!=0){
        mpz_clear(key);
        mpz_clear(n);
        mpz_clear(msg);
        return EXIT_FAILURE;
    }
    
    mpz_t r;
    mpz_init(r);
    mpz_powm(r,msg,key,n);
    
    gmp_printf("%#Zx\n",r);
    
    mpz_clear(key);
    mpz_clear(n);
    mpz_clear(msg);
    mpz_clear(r);
    
    return EXIT_SUCCESS;
}

/*
	Triviální faktorizace. Pokud je ůspěšná, faktor je vypsán na stdout.
	Vrací EXIT_SUCCESS nebo EXIT_FAILURE.
	n číslo pro faktorizaci
*/
int trivial_break(mpz_t n){
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_t i;
    mpz_init(i);
    
    for(mpz_set_ui(i,2ul);mpz_cmp_ui(i,1000000ul)<0;mpz_add_ui(i,i,1ul)){
        mpz_mod(tmp,n,i);
        if(mpz_cmp_ui(tmp,0ul)==0){
            gmp_printf("%#Zx\n",i);
            mpz_clear(tmp);
            mpz_clear(i);
            return EXIT_SUCCESS;
        }
    }
    
    mpz_clear(tmp);
    mpz_clear(i);
    return EXIT_FAILURE;
}

/*
	Nalezne minimum ze dvou čísel.
	r výsledné minimum mezi a, b
	a číslo
	b číslo
*/
void min(mpz_t r, mpz_t a, mpz_t b){
	if(mpz_cmp(a,b)<0){
		mpz_set(r,a);
	}
	else{
		mpz_set(r,b);
	}
}

/*
	Faktorizace pomocí brentova algoritmu. Pokud je ůspěšná, faktor je vypsán
	na stdout.
	Vrací EXIT_SUCCESS nebo EXIT_FAILURE.
	n číslo pro faktorizaci
*/
int brent_break(mpz_t n){
	if(mpz_even_p(n)){
		printf("%#x\n",2);
		return EXIT_SUCCESS;
	}
	
	gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,seed());
    
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_t top;
    mpz_init(top);
    mpz_sub_ui(top,n,1ul);
    
    mpz_t x;
    mpz_init(x);
    
    mpz_t y;
    mpz_init(y);
    mpz_urandomm(y,state,top);
    mpz_add_ui(y,y,1ul);
    
    mpz_t ys;
    mpz_init(ys);
    
    mpz_t c;
    mpz_init(c);
    mpz_urandomm(c,state,top);
    mpz_add_ui(c,c,1ul);
    
    mpz_t m;
    mpz_init(m);
    mpz_urandomm(m,state,top);
    mpz_add_ui(m,m,1ul);
    
    mpz_t g;
    mpz_init(g);
    mpz_set_ui(g,1ul);
    
    mpz_t r;
    mpz_init(r);
    mpz_set_ui(r,1ul);
    
    mpz_t q;
    mpz_init(q);
    mpz_set_ui(q,1ul);
    
    mpz_t i;
    mpz_init(i);
    
    mpz_t k;
    mpz_init(k);
    
    while(mpz_cmp_ui(g,1ul)==0){
    	mpz_set(x,y);
    	for(mpz_set_ui(i,0ul);mpz_cmp(i,r)<0;mpz_add_ui(i,i,1ul)){
    		mpz_powm_ui(y,y,2ul,n);
    		mpz_add(y,y,c);
    		mpz_mod(y,y,n);
    	}
    	mpz_set_ui(k,0ul);
    	while(mpz_cmp(k,r)<0 && mpz_cmp_ui(g,1ul)==0){
    		mpz_set(ys,y);
    		
    		mpz_sub(top,r,k);
    		min(top,m,top);
    		for(mpz_set_ui(i,0ul);mpz_cmp(i,top)<0;mpz_add_ui(i,i,1ul)){
    			mpz_powm_ui(y,y,2ul,n);
				mpz_add(y,y,c);
				mpz_mod(y,y,n);
				
				mpz_sub(tmp,x,y);
				mpz_abs(tmp,tmp);
				mpz_mul(q,q,tmp);
				mpz_mod(q,q,n);
    		}
    		gcd(g,q,n);
    		mpz_add(k,k,m);
    	}
    	mpz_mul_ui(r,r,2ul);
    }
    
    if(mpz_cmp(g,n)==0){
    	for(;;){
    		mpz_mul(ys,ys,ys);
    		mpz_mod(ys,ys,n);
    		mpz_add(ys,ys,c);
    		mpz_mod(ys,ys,n);
    		
    		mpz_sub(tmp,x,ys);
    		mpz_abs(tmp,tmp);
    		gcd(g,tmp,n);
    		
    		if(mpz_cmp_ui(g,1ul)>0){
    			break;
    		}
    	}
    }
    
    gmp_printf("%#Zx\n",g);
    
    mpz_clear(top);
    mpz_clear(tmp);
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(ys);
    mpz_clear(c);
    mpz_clear(m);
    mpz_clear(g);
    mpz_clear(r);
    mpz_clear(q);
    mpz_clear(i);
    mpz_clear(k);
    gmp_randclear(state);
	
	return EXIT_SUCCESS;
}

/*
	Faktorizace pomocí pollardova p-1 algoritmu. Pokud je ůspěšná, faktor je vypsán
	na stdout.
	Vrací EXIT_SUCCESS nebo EXIT_FAILURE.
	n číslo pro faktorizaci
*/
/*
int pollardp1_break(mpz_t n){
    mpz_t a;
    mpz_init(a);
    
    mpz_t r;
    mpz_init(r);
    
    mpz_t a_to_r_fac;
    mpz_init(a_to_r_fac);
    
    mpz_t d;
    mpz_init(d);
    
    mpz_t tmp;
    mpz_init(tmp);
    gmp_printf("n=%#Zx\n",n);
    
    for(mpz_set_ui(a,2ul);;mpz_add_ui(a,a,1ul)){
		gmp_printf("a=%#Zx\n",a);
    	
    	gcd(tmp,a,n);
    	if(mpz_cmp_ui(tmp,1ul)>0){
    		gmp_printf("%#Zx\n",tmp);
    		return EXIT_SUCCESS;
    	}
    	
    	mpz_set(a_to_r_fac,a);
    	for(mpz_set_ui(r,2ul);;mpz_add_ui(r,r,1ul)){
    		//gmp_printf("r=%#Zx\n",r);
    		mpz_powm(a_to_r_fac,a_to_r_fac,r,n);
    		mpz_sub_ui(tmp,a_to_r_fac,1ul);
    		gcd(d,tmp,n);
    		//gmp_printf("t=%#Zx\n",tmp);
    		if(mpz_cmp(d,n)==0){
    			break;
    		}
    		else if(mpz_cmp_ui(d,1ul)>0){
    			gmp_printf("%#Zx\n",d);
    			return EXIT_SUCCESS;
    			
    		}
    	}	
    }
    
    mpz_clear(a);
    mpz_clear(r);
    mpz_clear(a_to_r_fac);
    mpz_clear(d);
    mpz_clear(tmp);
    
    return EXIT_FAILURE;
}
*/

/*
	Faktorizace pomocí fermatova algoritmu. Pokud je ůspěšná, faktor je vypsán
	na stdout.
	Vrací EXIT_SUCCESS nebo EXIT_FAILURE.
	n číslo pro faktorizaci
*/
/*
int fermat_break(mpz_t n){
    mpz_t a;
    mpz_init(a);
    
    mpz_sqrt(a,n);
    mpz_add_ui(a,a,1ul);
    
    mpz_t b2;
    mpz_init(b2);
    mpz_mul(b2,a,a);
    mpz_sub(b2,b2,n);
    
    while(!mpz_perfect_square_p(b2)){
        mpz_add_ui(a,a,1ul);
        mpz_mul(b2,a,a);
        mpz_sub(b2,b2,n);
    }
    
    mpz_sqrt(b2,b2);
    mpz_sub(a,a,b2);
    gmp_printf("%#Zx\n",a);
    
    mpz_clear(a);
    mpz_clear(b2);
    return EXIT_SUCCESS;
}
*/

/*
	Faktorizace čísla n. Pokud je ůspěšná, faktor je vypsán
	na stdout.
	Vrací EXIT_SUCCESS nebo EXIT_FAILURE.
	n číslo pro faktorizaci
*/
int break_cypher(Args args){
    int r=EXIT_FAILURE;
    int ctrl=1;
    
    //číslo pro faktorizaci
    mpz_t n;
    mpz_init(n);
    ctrl=mpz_set_str(n,args.n,0);
    if(ctrl!=0){
        mpz_clear(n);
        return EXIT_FAILURE;
    }
    
    //triviální faktorizace
    r=trivial_break(n);
    
    //brentův algoritmus
    if(r==EXIT_FAILURE){
        r=brent_break(n);
    }
    
    mpz_clear(n);
    
    return r;
}

/*
	Hlavní funkce.
	argc počet argumentů
	argv argumenty
*/
int main(int argc, char** argv){
    Args args=get_args(argc, argv);
    int r=EXIT_FAILURE;
    
    switch(args.a){
        case ERROR:r=EXIT_FAILURE;break;
        case GENERATE:r=generate(args);break;
        case ENCRYPT:r=encrypt(args);break;
        case DECRYPT:r=encrypt(args);break;
        case BREAK:r=break_cypher(args);break;
    }
    
    return r;
}

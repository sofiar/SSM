# prueba super chiquita para ver si llamar rexp() desde Cpp demora
library(tidyverse); library(microbenchmark)

# Funcion que sortea n valores exp(1) usando rexp de R
cppFunction('
            NumericVector fexp(double n) {
             	Function rexp("rexp");
              NumericVector x(n);
               x = rexp(n, 1.0);
 		          return x;
            }'
            )

# Funcion que sortea n valores exp(1) usando rexp de Cpp
cppFunction('
            NumericVector fexp_cpp(double n) {
              NumericVector x(n);
               x = rexp(n, 1.0);
 		          return x;
            }'
)

# Miramos los tiempos: 
m <- microbenchmark( list=list( f1_Ryc = {fexp(1e4)},  # llamar R desde Cpp
                                f2_Cpp={fexp_cpp(1e4)},  # solo Cpp
                                f3_R={rexp(1e4, 1)} )  # solo R
                     )

m %>% ggplot(aes(x=expr, y=time)) + geom_boxplot() + 
  scale_y_log10() + 
  labs(title = 'no cambia un corno')

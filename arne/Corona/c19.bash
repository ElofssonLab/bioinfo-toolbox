#!/bin/bash -x

regions=Stockholm Östergötland "Västra Götaland" Skåne Sörmland Uppsala Jönköping Örebro Dalarna Halland Gävleborg Västmanland Jämtland Norrbotten Västerbotten Västernorrland Värmland Kronoberg Kalmar Blekinge Gotland Okänt

dir=~/Desktop/Corona/c19/data/

for i in Stockholm Östergötland "Västra Götaland" Skåne Sörmland Uppsala Jönköping Örebro Dalarna Halland Gävleborg Västmanland Jämtland Norrbotten Västerbotten Västernorrland Värmland Kronoberg Kalmar Blekinge Gotland Okänt
do
    curl https://c19.se/Sweden/$i | grep -E "name:|data:|categories:" > $dir/$i.txt
done
# datum, Fall, Fall idag, döda, Total IVA


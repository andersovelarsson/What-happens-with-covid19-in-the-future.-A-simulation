# Vad händer med Covid19 i framtiden. 


## S - I - R - D model

Susceptibles   S(t)

Infected       I(t)

Recovered      R(t)  

Death          D(t)

Population N(t) = S(t) +  I(t) + R(t) - D(t)

Initial condition
N(0)  =  No  = 10 000 001

S(0)  =  So  = 10 000 000

I(0)  =  Io  = 1

R(0)  =  Ro  = 0

dS/dt = -k01 * I * S                       (0)

dI/dt =  k01 * I * S - k12 * I - k13 * I   (1)

dR/dt =  k12 * I                           (2)

dD/dt =  k13 * I                           (3)

Det den matamatiska modellen säger är att Covid-19 pandemin troligen handlar om dig och mig och hur uthålliga vi kan vara. Frågan är om det kommer en andra och tredje våg av epedemin. Modellen kan inte svara exakt på denna fråga, men modellen visar kopplingen mellan vårat betende och utfall. Svaret på frågan verkar bero på om vi går tillbaka och lever som vi gjorde innan (Ro = 2.5) pandemin då kommer andra vågen troligen bli mycket starakre än första vågen, kanske 2-3 gånger starkare. Forsätter vi leva som vi gör nu (början av juni 2020, Ro < 1.4) då kommer troligen ingen andra våg, men vi kommer att behöva vara uthålliga. Troligen kommer vi behöva hålla ut till vi har ett vacsin. 

![](/Doc/SimulationResult.png)


![](/Doc/SimulationResult2.png)

En relativ enkel modell för hur Covid19 sprids i samhället. Modellen kan förstås inte förutsäga hur spridningen av Covid19 kommer att se ut i mista detalj. Det finns några saker som är svårt att förutsäga och som spelar en stor roll för resultatet vad som händer i framtiden med Covid19. Det första är att det är svårt att förutse hur vi mäniskor kommer att ändra sitt betende när det gäller smittspridning. Det andra är viruset i sig själv kan öka eller minska i sin förmåga att sprida sig med tiden. Dessa två fakroer påverkar det så kallade reproduktionstalet Ro. Modellen innehåller också vissa antagande för att förenkla modellen. Modellkomplexiteten är en kompromiss mellan något som är lätt att förstå, hur väl modellen fångar verkligheten samt tillgång till verklig data. Modellen fokuserar på de viktigaste fenomen i modellen snarare än alla detljer.
Även fast man inte kan estimera, med stor nogrannhet, vad som kommer att hända så är denna simulering ett verktyg för att titta på olika scenarier för utvecklingen av Covid19 i samhället.

Modellen är anpassad utfrån data från FHM

## License
[MIT](https://choosealicense.com/licenses/mit/)



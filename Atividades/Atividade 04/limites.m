function vet = limites( vet, Min, Max)
  vet(vet < Min) = Min;
  vet(vet > Max) = Max;
end
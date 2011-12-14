function cost=costOne(len,gc)
   %attribute	max	mean	min
   %gc   	.45	.23	.13
   %length	30864	 25573	112
   %
   %lenght is good
   %high gc is good
   %
   cost=(len/25573  + gc/.23).^-1;

   %histogram centered near cost of 100
end


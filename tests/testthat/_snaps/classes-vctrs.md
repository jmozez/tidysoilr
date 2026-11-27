# soil_profile() builds a classed, sorted tibble

    Code
      print(prof, n = 3, width = 80)
    Output
      # A <soil_profile>: 150 x 25
      # Profiles:         32 (id: profile)
      # Depths:           0 to 402 (top, bottom)
        profile dataset licence        wosis_id  year   lon   lat wrb_group usda_order
        <chr>   <chr>   <chr>             <int> <int> <dbl> <dbl> <chr>     <chr>     
      1 82P0793 US-NCSS US Public Dom~  1052498  1982  145. -6.03 <NA>      <NA>      
      2 82P0793 US-NCSS US Public Dom~  1052498  1982  145. -6.03 <NA>      <NA>      
      3 82P0793 US-NCSS US Public Dom~  1052498  1982  145. -6.03 <NA>      <NA>      
      # i 147 more rows
      # i 16 more variables: horizon <chr>, top <dbl>, bottom <dbl>, sand <dbl>,
      #   silt <dbl>, clay <dbl>, coarse_frag <dbl>, oc <dbl>, total_n <dbl>,
      #   bulk_density <dbl>, ph_water <dbl>, ph_kcl <dbl>, cec <dbl>, ecec <dbl>,
      #   theta_33 <dbl>, theta_1500 <dbl>

# soil_spectra() reads wavelengths and validates

    Code
      print(spec, n = 2, width = 80)
    Output
      # A <soil_spectra>: 120 x 112
      # Spectra:          105 bands, 400 to 2480 (reflectance)
      # Preprocessing:    none
        sample profile depth   lon   lat  clay    oc nm_400 nm_420 nm_440 nm_460
        <chr>  <chr>   <chr> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
      1 N001   P01     0-15   149. -4.62  29.9  12    0.132  0.139  0.148  0.156
      2 N002   P01     15-30  149. -4.62  26.7  23.6  0.109  0.117  0.124  0.132
      # i 118 more rows
      # i 101 more variables: nm_480 <dbl>, nm_500 <dbl>, nm_520 <dbl>, nm_540 <dbl>,
      #   nm_560 <dbl>, nm_580 <dbl>, nm_600 <dbl>, nm_620 <dbl>, nm_640 <dbl>,
      #   nm_660 <dbl>, nm_680 <dbl>, nm_700 <dbl>, nm_720 <dbl>, nm_740 <dbl>,
      #   nm_760 <dbl>, nm_780 <dbl>, nm_800 <dbl>, nm_820 <dbl>, nm_840 <dbl>,
      #   nm_860 <dbl>, nm_880 <dbl>, nm_900 <dbl>, nm_920 <dbl>, nm_940 <dbl>,
      #   nm_960 <dbl>, nm_980 <dbl>, nm_1000 <dbl>, nm_1020 <dbl>, ...


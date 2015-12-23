function gmt = calculate_GMT(date)

day = mod(date, 1);
gmt = day*24;


The program `ascii2bin` converts an ASCII file in the format of EMmod into a binary EMmod file.

Compilation:

```
gcc ascii2bin.c -DDOUBLE -o ascii2bin
```

The flag `-DDOUBLE` ensures double precision is used when writing the binary file. 

Usage:

```
 ./ascii2bin file_in number_of_elements file_out
```

where
- `file_in` is the name of the ASCII input file.
- `number_of_elements` is an integer specifying how many elements need to be read. If the dataset contains for example 512 times 512 datapoints, `number_of_elements` is 262144.
- `file_out` is the name of the binary output file.

module Not_gate
(
 a,
 b,
 out
 );

   input a;
   input b;

   output out;

   assign out = a not b;
 
endmodule // not_gate


module nor_gate
(
 a,
 b,
 out
 );

   input a;
   input b;

   output out;

   assign out = a not or b;
 
endmodule // nor_gate

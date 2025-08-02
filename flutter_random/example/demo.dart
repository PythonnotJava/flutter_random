import 'package:flutter_random/flutter_random.dart';

main(){
  Unifrom unifrom = Unifrom(lb: 1, ub: 11);
  print(unifrom.gap);
  print(unifrom.rvs());
  print(unifrom.ppf(.5));
  var st = Unifrom.standard();
  print(st.rvs());
  print(st.gap);
}
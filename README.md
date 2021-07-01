# LDPC
special topic LDPC<br/>
implement algorithms including log SPA, min sum, min sum-C<br/>
<br/>
- LDPC_SPA1 uses the variables Q, q, and r. This method first came up of my mind when I tried to programming SPA.<br/>
- LDPC_SPA2 uses the variables Q, r, and prev_r. The implementation involves the boxsum operator, and this helps LDPC_SPA2 has lower time complexity than LDPC_SPA1. <br/>
- Min sum-C is the min-sum algorithm with constant correction factor. I chose c=0.1 since it achieved better performance at c=0.1, but I don't know the reason.
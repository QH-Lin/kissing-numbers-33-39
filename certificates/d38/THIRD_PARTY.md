# Third-party material in the 38-dimensional certificate

The 38-dimensional configuration adds new points to a published one. Only `equator_lines_6855.txt` is new.

| file | source | license |
|---|---|---|
| `base-data/*` (Golay basis, class labels of the 98280 minimal Leech lines, 14-dimensional directions and triples, axis permutation, rational rotation) | Alexey Kravatskiy, [KissingNumbers](https://github.com/alexlegeartis/KissingNumbers/tree/cf14c5ef4db6059bbf2e570e3f0ce24edfda243e/verifications/improved/dim38-leech-large-codimension), revision cf14c5e, the 591,612-point configuration; packaged as `constructions/d38/base-data/` in the [Qiushi Engine repository](https://github.com/Oxelra-AI/Qiushi-Engine-Kissing-Number-Research) (revision bbb7aaf), from which the files are copied byte for byte | MIT (Kravatskiy), notice below; the Qiushi packaging is CC BY 4.0 |
| `../../verify/d38_verify.py` | adapted from `tools/verify_dimension38.py` of the Qiushi Engine repository, revision bbb7aaf; the changes are listed in its docstring | MIT, notice below |

The idea of adding norm-48 Leech vectors on the empty equator of Kravatskiy's configuration, and the 144-line version of it (591,900 points), are due to S. Yang et al., *Large-scale autonomous discovery of kissing number constructions*, [arXiv:2609.35051](https://arxiv.org/abs/2609.35051), §7.

## KissingNumbers (MIT License)

```
MIT License

Copyright (c) 2026 Alexey Kravatskiy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Qiushi Engine Kissing Number Research (software: MIT License)

```
MIT License

Copyright (c) 2026 Qiushi Engine contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

The Qiushi Engine repository offers its original non-software data under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/); the dimension-39 dense layer files `certificates/d39/dense_*.txt` come from there as well.

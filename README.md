fzero
=====

Find a zero of a univariate function.  This is a JavaScript port of [Jaroslav Hajek](highegg@gmail.com)'s [fzero implementation](https://fossies.org/dox/octave-4.0.0/fzero_8m_source.html) in [Octave](https://www.gnu.org/software/octave/).

This is essentially the ACM algorithm 748: Enclosing Zeros of Continuous Functions due to Alefeld, Potra and Shi, ACM Transactions on Mathematical Software, Vol. 21, No. 3, September 1995.  Although the workflow should be the same, the structure of the algorithm has been transformed non-trivially; instead of the authors' approach of sequentially calling building blocks subprograms we implement here a FSM version using one interior point determination and one bracketing per iteration, thus reducing the number of temporary variables and simplifying the algorithm structure.  Further, this approach reduces the need for external functions and error handling.  The algorithm has also been slightly modified.

Usage
-----
```
$ npm install fzero
```
To use fzero in Node.js, simply require it:
```javascript
var fzero = require("fzero");
```
A minified, browserified file `dist/fzero.min.js` is included for use in the browser.  Including this file attaches a `fzero` object to `window`:
```html
<script src="dist/fzero.min.js" type="text/javascript"></script>
```

Tests
-----
Unit tests are included in `test/`, and can be run using npm:
```
$ npm test
```

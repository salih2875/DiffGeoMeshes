let gauss = [];
let mean = [];
let firstcurv = [];
let secondcurv = [];
let gcolors = [];
let mcolors = [];
let fcolors = [];
let scolors = [];

Module.onRuntimeInitialized = async function () {
  document
    .getElementById("compute")
    .addEventListener("click", async function () {
      let timea = performance.now();

      document.getElementById("status").textContent = "...computing";

      let meshfile = document.getElementById("url").value;
      console.log(meshfile);
      const response = await fetch(meshfile);
      const body = await response.text();
      let flag = checkformat(body);
      let vertices = [];
      let indexes = [];
      if (flag == true) {
        parseoff(body, vertices, indexes);
      } else {
        parseobj(body, vertices, indexes);
      }

      let v = new Float64Array(vertices);
      let t = new Int32Array(indexes);

      console.log("wasm running");
      var len = v.length;
      var bytes_per_element = v.BYTES_PER_ELEMENT;
      var offset = _malloc(len * bytes_per_element);

      var len1 = t.length;
      var bytes_per_element1 = t.BYTES_PER_ELEMENT;
      var offset2 = _malloc(len1 * bytes_per_element1);

      let gauss_ptr = Module._malloc((len / 3) * bytes_per_element);
      let mean_ptr = Module._malloc((len / 3) * bytes_per_element);
      let first_ptr = Module._malloc((len / 3) * bytes_per_element);
      let second_ptr = Module._malloc((len / 3) * bytes_per_element);

      let gcolor = Module._malloc(len * bytes_per_element);
      let mcolor = Module._malloc(len * bytes_per_element);
      let fcolor = Module._malloc(len * bytes_per_element);
      let scolor = Module._malloc(len * bytes_per_element);

      Module.HEAPF64.set(v, offset / bytes_per_element);
      Module.HEAP32.set(t, offset2 / bytes_per_element1);

      _printgauss(
        len / 3,
        len1 / 3,
        offset,
        offset2,
        gauss_ptr,
        mean_ptr,
        first_ptr,
        second_ptr,
        gcolor,
        mcolor,
        fcolor,
        scolor
      );

      gauss = new Float64Array(Module.HEAPF64.buffer, gauss_ptr, len / 3);
      mean = new Float64Array(Module.HEAPF64.buffer, mean_ptr, len / 3);
      firstcurv = new Float64Array(Module.HEAPF64.buffer, first_ptr, len / 3);
      secondcurv = new Float64Array(Module.HEAPF64.buffer, second_ptr, len / 3);

      gcolors = new Float64Array(Module.HEAPF64.buffer, gcolor, len);
      mcolors = new Float64Array(Module.HEAPF64.buffer, mcolor, len);
      fcolors = new Float64Array(Module.HEAPF64.buffer, fcolor, len);
      scolors = new Float64Array(Module.HEAPF64.buffer, scolor, len);

      console.log(gauss);
      console.log(mean);
      console.log(firstcurv);
      console.log(secondcurv);

      _free(offset);
      _free(offset2);
      _free(gauss_ptr);
      _free(mean_ptr);
      _free(first_ptr);
      _free(second_ptr);

      document.getElementById("status").textContent = "Done!";

      let timeb = performance.now();

      let perf = timeb - timea;
      console.log(perf);
      // window.alert("Calculated!!\nSet scale and click render");
    });
};

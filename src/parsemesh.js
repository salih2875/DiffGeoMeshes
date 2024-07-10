let vertices = [];
let indexes = [];

function checkformat(data) {
  let pattern = /OFF/g;
  let splitData = data.split("\n");
  for (var i = 0, l = splitData.length; i < l; i++) {
    if (pattern.exec(splitData[i]) != null) {
      return true;
    }
  }

  return false;
}


// https://gist.github.com/phoboslab/5964733
function parseobj(data, vertices, indexes) {
  let pattern, result;

  // v float float float
  pattern = /v( +[\d|\.|\+|\-|e]+)( [\d|\.|\+|\-|e]+)( [\d|\.|\+|\-|e]+)/g;

  while ((result = pattern.exec(data)) != null) {
    // ["v 1.0 2.0 3.0", "1.0", "2.0", "3.0"]

    vertices.push(
      parseFloat(result[1]),
      parseFloat(result[2]),
      parseFloat(result[3])
    );
  }

  var splitData = data.split("\no ");

  for (var i = 0, l = splitData.length; i < l; i++) {
    let object = splitData[i];

    // f vertex vertex vertex ...
    pattern = /f( +[\d]+)( [\d]+)( [\d]+)( [\d]+)?/g;
    while ((result = pattern.exec(object)) != null) {
      // ["f 1 2 3", "1", "2", "3", undefined]
      if (result[4] === undefined) {
        indexes.push(
          parseInt(result[1]),
          parseInt(result[2]),
          parseInt(result[3])
        );
      } else {
        indexes.push(
          parseInt(result[1]),
          parseInt(result[2]),
          parseInt(result[3]),
          parseInt(result[4])
        );
      }
    }

    // f vertex/uv vertex/uv vertex/uv ...
    pattern =
      /f( +([\d]+)\/([\d]+))( ([\d]+)\/([\d]+))( ([\d]+)\/([\d]+))( ([\d]+)\/([\d]+))?/g;

    while ((result = pattern.exec(object)) != null) {
      // ["f 1/1 2/2 3/3", " 1/1", "1", "1", " 2/2", "2", "2", " 3/3", "3", "3", undefined, undefined, undefined]

      if (result[10] === undefined) {
        indexes.push(
          parseInt(result[2]) - 1,
          parseInt(result[5]) - 1,
          parseInt(result[8]) - 1
        );
      } else {
        indexes.push(
          parseInt(result[2]) - 1,
          parseInt(result[5]) - 1,
          parseInt(result[8]) - 1,
          parseInt(result[11]) - 1
        );
      }
    }

    // f vertex/uv/normal vertex/uv/normal vertex/uv/normal ...

    pattern =
      /f( +([\d]+)\/([\d]+)\/([\d]+))( ([\d]+)\/([\d]+)\/([\d]+))( ([\d]+)\/([\d]+)\/([\d]+))( ([\d]+)\/([\d]+)\/([\d]+))?/g;

    while ((result = pattern.exec(object)) != null) {
      // ["f 1/1/1 2/2/2 3/3/3", " 1/1/1", "1", "1", "1", " 2/2/2", "2", "2", "2", " 3/3/3", "3", "3", "3", undefined, undefined, undefined, undefined]

      if (result[13] === undefined) {
        indexes.push(
          parseInt(result[2]) - 1,
          parseInt(result[6]) - 1,
          parseInt(result[10]) - 1
        );
      } else {
        indexes.push(
          parseInt(result[2]) - 1,
          parseInt(result[6]) - 1,
          parseInt(result[10]) - 1,
          parseInt(result[14]) - 1
        );
      }
    }

    // f vertex//normal vertex//normal vertex//normal ...

    pattern =
      /f( +([\d]+)\/\/([\d]+))( ([\d]+)\/\/([\d]+))( ([\d]+)\/\/([\d]+))( ([\d]+)\/\/([\d]+))?/g;

    while ((result = pattern.exec(object)) != null) {
      // ["f 1//1 2//2 3//3", " 1//1", "1", "1", " 2//2", "2", "2", " 3//3", "3", "3", undefined, undefined, undefined]

      if (result[10] === undefined) {
        indexes.push(
          parseInt(result[2]) - 1,
          parseInt(result[5]) - 1,
          parseInt(result[8]) - 1
        );
      } else {
        indexes.push(
          parseInt(result[2]) - 1,
          parseInt(result[5]) - 1,
          parseInt(result[8]) - 1,
          parseInt(result[11]) - 1
        );
      }
    }
  }
}

function parseoff(data, vertices, indexes) {
  let pattern, result;

  var splitData = data.split("\n");
  // let sizes = splitData[1];
  // console.log(sizes);
  let nv, nf;

  pattern = /([\d]+)( [\d]+)( [\d]+)/g;

  while ((result = pattern.exec(splitData[1])) != null) {
    nv = parseInt(result[1]);
    // console.log(nv);
    nf = parseInt(result[2]);
    // console.log(nf);
  }

  // parse floats
  for (var i = 2, l = nv + 2; i < l; i++) {
    if (i == nv - 1) {
      console.log("parsing is done");
    }
    let line = splitData[i];
    // console.log(line);
    pattern = /([\d|\.|\+|\-|e]+)( [\d|\.|\+|\-|e]+)( [\d|\.|\+|\-|e]+)/g;
    while ((result = pattern.exec(line)) != null) {
      vertices.push(
        parseFloat(result[1]),
        parseFloat(result[2]),
        parseFloat(result[3])
      );
    }
  }

  // parse indexes
  for (var i = nv + 2, l = nf + nv + 2; i < l; i++) {
    // console.log(nv+2);
    let object = splitData[i];

    // console.log(object);

    // f vertex vertex vertex ...
    pattern = /3( +[\d]+)( [\d]+)( [\d]+)/g;
    while ((result = pattern.exec(object)) != null) {
      // ["f 1 2 3", "1", "2", "3", undefined]
      indexes.push(
        parseInt(result[1]),
        parseInt(result[2]),
        parseInt(result[3])
      );
    }
  }
}

async function load(url = "") {
  const response = await fetch(url);
  const body = await response.text();

  let flag = checkformat(body);

  let vertices = [];
  let indexes = [];
  if (flag == true) {
    parseoff(body, vertices, indexes);
  } else {
    parseobj(body, vertices, indexes);
  }
}

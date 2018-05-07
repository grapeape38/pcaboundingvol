function jacobi(m) {
  const MAX_SWEEPS = 1000;
  const epsilon = 0.00001;
  var m11 = m[0][0];
  var m12 = m[0][1];
  var m13 = m[0][2];
  var m22 = m[1][1];
  var m23 = m[1][2];
  var m33 = m[2][2];
  r = [[1,0,0],[0,1,0],[0,0,1]];
  for (var a = 0; a < MAX_SWEEPS; a++) {
    // Exit if off-diagonal entries small enough.
    if ((Math.abs(m12) < epsilon) && (Math.abs(m13) < epsilon) &&
      (Math.abs(m23) < epsilon)) break;
    // Annihilate (1,2) entry.
    if (m12 != 0.0) {
      var u = (m22 - m11) * 0.5 / m12;
      var u2 = u * u;
      var u2p1 = u2 + 1.0;
      var t = (u2p1 != u2) 
        ? ((u < 0.0) ? -1.0 : 1.0) * (Math.sqrt(u2p1) - Math.abs(u)) 
        : 0.5 / u;
      var c = 1.0 / Math.sqrt(t * t + 1.0);
      var s = c * t;
      m11 -= t * m12;
      m22 += t * m12;
      m12 = 0.0;
      var temp = c * m13 - s * m23;
      m23 = s * m13 + c * m23;
      m13 = temp;
      for (var i = 0; i < 3; i++) {
        var temp = c * r[i][0] - s * r[i][1];
        r[i][1] = s * r[i][0] + c * r[i][1];
        r[i][0] = temp;
      }
    }
    // Annihilate (1,3) entry.
    if (m13 != 0.0) {
      var u = (m33 - m11) * 0.5 / m13;
      var u2 = u * u;
      var u2p1 = u2 + 1.0;
      var t = (u2p1 != u2) ?
      ((u < 0.0) ? -1.0 : 1.0) * (Math.sqrt(u2p1) - Math.abs(u))
      : 0.5 / u;
      var c = 1.0 / Math.sqrt(t * t + 1.0);
      var s = c * t;
      m11 -= t * m13;
      m33 += t * m13;
      m13 = 0.0;
      var temp = c * m12 - s * m23;
      m23 = s * m12 + c * m23;
      m12 = temp;
      for (var i = 0; i < 3; i++) {
        var temp = c * r[i][0] - s * r[i][2];
        r[i][2] = s * r[i][0] + c * r[i][2];
        r[i][0] = temp;
      }
    }
    // Annihilate (2,3) entry.
    if (m23 != 0.0) {
      var u = (m33 - m22) * 0.5 / m23;
      var u2 = u * u;
      var u2p1 = u2 + 1.0;
      var t = (u2p1 != u2) ?
      ((u < 0.0) ? -1.0 : 1.0) * (Math.sqrt(u2p1) - Math.abs(u))
      : 0.5 / u;
      var c = 1.0 / Math.sqrt(t * t + 1.0);
      var s = c * t;
      m22 -= t * m23;
      m33 += t * m23;
      m23 = 0.0;
      var temp = c * m12 - s * m13;
      m13 = s * m12 + c * m13;
      m12 = temp;
      for (var i = 0; i < 3; i++) {
        var temp = c * r[i][1] - s * r[i][2];
        r[i][2] = s * r[i][1] + c * r[i][2];
        r[i][1] = temp;
      }
    }
  }
  var lambda = [m11,m22,m33];
  return {'lambda': lambda, 'r': r};
}

var min = -6;
var max = 6;
function randDouble() {
  return Math.random() * (max-min) + min;
}

function getMeanVector(rows) {
  var len = rows.length;
  var m = [0,0,0]
  for (var i = 0; i < len; i++) {
    for (var j = 0; j < 3; j++) {
      m[j] += rows[i][j];
    }
  }
  for (var i = 0; i < 3; i++) m[i] /= len;
  return m;
}

function getCovarianceMatrix(P, m) {
  C = [[0,0,0],[0,0,0],[0,0,0]];
  for (var a = 0; a < P.length; a++) {
    var vec = P[a];
    for (var i = 0; i < 3; i++) {
      for(var j = 0; j < 3; j++) {
        C[i][j] += (vec[i] - m[i]) * (vec[j] - m[j]);
      }
    }
  }
  for (var i = 0; i < 3; i++)
    for (var j = 0; j < 3; j++)
      C[i][j] /= P.length; 
  return C;
}

function dot(u, v) {
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; 
}

function find3x3Inverse(A) {
  var a = A[0][0], b = A[0][1], c = A[0][2],
      d = A[1][0], e = A[1][1], f = A[1][2],
      g = A[2][0], h = A[2][1], i = A[2][2];
  var det = a*(e*i-f*h) - b*(d*i - f*g) + c*(d*h - e*g);
  var s = 1 / det;
  var inv = [[(e*i - f*h) * s, (c*h - b*i) * s, (b*f - c*e)*s],
             [(f*g - d*i) * s, (a*i - c*g) * s, (c*d - a*f)*s],
             [(d*h - e*g) * s, (b*g - a*h) * s, (a*e - b*d)*s]];
  return inv;
}

function findCubeVertex(p1, p2, p3, d1, d2, d3) {
  var A = [[p1[0],p1[1],p1[2]],
           [p2[0],p2[1],p2[2]],
           [p3[0],p3[1],p3[2]]];
  var inv = find3x3Inverse(A);
  var vertex = [inv[0][0] * -d1 + inv[0][1] * -d2 + inv[0][2] * -d3,
                inv[1][0] * -d1 + inv[1][1] * -d2 + inv[1][2] * -d3,
                inv[2][0] * -d1 + inv[2][1] * -d2 + inv[2][2] * -d3];
  return vertex;
}

function findVertices(eigenMatrix, extrema) {
  var R = [eigenMatrix[0][0], eigenMatrix[1][0], eigenMatrix[2][0]];
  var S = [eigenMatrix[0][1], eigenMatrix[1][1], eigenMatrix[2][1]];
  var T = [eigenMatrix[0][2], eigenMatrix[1][2], eigenMatrix[2][2]];
  var Rneg = [-eigenMatrix[0][0], -eigenMatrix[1][0], -eigenMatrix[2][0]];
  var Sneg = [-eigenMatrix[0][1], -eigenMatrix[1][1], -eigenMatrix[2][1]];
  var Tneg = [-eigenMatrix[0][2], -eigenMatrix[1][2], -eigenMatrix[2][2]];
  var v1 = findCubeVertex(R, S, T, -extrema[0], -extrema[1], -extrema[2]);
  var v2 = findCubeVertex(R, S, Tneg, -extrema[0], -extrema[1], extrema[5]);
  var v3 = findCubeVertex(R, Sneg, T, -extrema[0], extrema[4], -extrema[2]);
  var v4 = findCubeVertex(R, Sneg, Tneg, -extrema[0], extrema[4], extrema[5]);
  var v5 = findCubeVertex(Rneg, S, T, extrema[3], -extrema[1], -extrema[2]);
  var v6 = findCubeVertex(Rneg, S, Tneg, extrema[3], -extrema[1], extrema[5]);
  var v7 = findCubeVertex(Rneg, Sneg, T, extrema[3], extrema[4], -extrema[2]);
  var v8 = findCubeVertex(Rneg, Sneg, Tneg, extrema[3], extrema[4], extrema[5]);
  return [v1,v2,v3,v4,v5,v6,v7,v8];
}

function findExtrema(P, eigenMatrix) {
  var R = [eigenMatrix[0][0], eigenMatrix[1][0], eigenMatrix[2][0]];
  var S = [eigenMatrix[0][1], eigenMatrix[1][1], eigenMatrix[2][1]];
  var T = [eigenMatrix[0][2], eigenMatrix[1][2], eigenMatrix[2][2]];
  var maxR = -100, maxS = -100, maxT = -100, minR = 100, minS = 100, minT = 100;
  for (var i = 0; i < P.length; i++) {
    if (dot(R, P[i]) > maxR)
      maxR = dot(R, P[i])
    if (dot(S, P[i]) > maxS)
      maxS = dot(S, P[i])
    if (dot(T, P[i]) > maxT)
      maxT = dot(T, P[i])
    if (dot(R, P[i]) < minR)
      minR = dot(R, P[i])
    if (dot(S, P[i]) < minS)
      minS = dot(S, P[i])
    if (dot(T, P[i]) < minT)
      minT = dot(T, P[i])
  }
  return [minR, minS, minT, maxR, maxS, maxT];
}

function unpack(rows, index) {
    return rows.map(function(row) 
    { return row[index]; });
}

function getCubeTrace(P, eigenMatrix) {
  var extrema = findExtrema(P, eigenMatrix);
  var vertices = findVertices(eigenMatrix, extrema);
  var cubeTrace = {
      type: "mesh3d",
      x: unpack(vertices, 0),
      y: unpack(vertices, 1), 
      z: unpack(vertices, 2), 
      i: [0, 1, 0, 1, 5, 5, 7, 7, 7, 7, 0, 2],
      j: [1, 2, 1, 4, 4, 6, 5, 1, 6, 2, 2, 4],
      k: [2, 3, 4, 5, 6, 7, 1, 3, 2, 3, 4, 6],
      opacity:0.0,
      color: 'blue',
    }
    return cubeTrace;
}

var data; 

function getTrace(rows) {
  return {
    x:unpack(rows, 0),  y: unpack(rows, 1), z: unpack(rows, 2), 
    mode: 'markers',
    marker: {
      size: 4,
      color: 'blue',
      opacity: 0.8
    },
    type: 'scatter3d',
    name: 'points'
  };
}

function genPoints() {
  var rows = []; 
  for (var i = 0; i < 100; i++) {
    x = randDouble();
    y = randDouble();
    z = randDouble();
    rows.push([x,y,z]);
  }
  return rows;
}

var layout = {
      dragmode: true,
      showlegend: false,
      margin: {
          l: 0,
          r: 0,
          b: 0,
          t: 0
       },
       scene: { 
          /*xaxis: {
            range: [-7, 7]
          },
          yaxis: {
            range: [-7, 7]
          },
          zaxis: {
            range: [-7, 7]
          }*/
      }
  };

function setLayout(cube_trace) {
  var minX = Math.min.apply(null, cube_trace.x);
  var minY = Math.min.apply(null, cube_trace.y);
  var minZ = Math.min.apply(null, cube_trace.z);
  var maxX = Math.max.apply(null, cube_trace.x);
  var maxY = Math.max.apply(null, cube_trace.y);
  var maxZ = Math.max.apply(null, cube_trace.z);
  var min = Math.min(minX,minY,minZ);
  var max = Math.max(maxX,maxY,maxZ);
  layout = {
      dragmode: true,
      showlegend: false,
      margin: {
          l: 0,
          r: 0,
          b: 0,
          t: 0
       },
       scene: { 
          xaxis: {
            range: [min, max]
          },
          yaxis: {
            range: [min, max]
          },
          zaxis: {
            range: [min, max]
          }
      }
  };
}

function loadData(rows) {
  var m = getMeanVector(rows);
  var C = getCovarianceMatrix(rows,m);
  var j = jacobi(C);
  var cube_trace = getCubeTrace(rows, j['r']);
  setLayout(cube_trace);
  var trace = getTrace(rows);
  data = [trace, cube_trace];
}

function plotPoints() {
  var choice = document.getElementById('plotSelect').value;
  var rows;
  if (choice == "random") {
    rows = genPoints();
  }
  else if (choice == "cat") {
    rows = cat;
  }
  else if (choice == "car") {
    rows = avent;
  }
  else {
    rows = house;
  }
  loadData(rows);
  Plotly.react('myDiv', data, layout);
  toggleBox();
}

loadData(genPoints());
Plotly.newPlot('myDiv', data, layout);

function toggleBox() {
  var showEigen = document.getElementById("showEigen").checked;
  var update = {
    opacity: showEigen ?  0.3 : 0.0
  };
  Plotly.restyle('myDiv', update, 1);
}

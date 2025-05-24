import { BufferAttribute, BufferGeometry, Color, CylinderGeometry, DirectionalLight, Euler, HemisphereLight, InstancedMesh, LineBasicMaterial, LineSegments, MeshLambertMaterial, Object3D, OrthographicCamera, PCFSoftShadowMap, Quaternion, Scene, Vector3, WebGLRenderer } from "https://esm.sh/three@0.169.0"


const { cos, min, PI, random: rnd, sin, sqrt } = Math
const TAU = PI * 2
const pr = min(2, devicePixelRatio)

const N = 1e3 * pr
const MAX_R = 0.05
const MIN_R = 0.006 / pr
const MAX_H = 0.1
const MIN_H = MIN_R * 2
const MAX_TESTS = N * 1e4


//////// SETUP

const scene = new Scene()
scene.background = new Color(0.01, 0.01, 0.01)

const camera = new OrthographicCamera()
camera.near = 0
camera.far = 4
set_camera_aspect(camera)
camera.position.set(0, 0, 2)

const renderer = new WebGLRenderer({ canvas })
renderer.setSize(innerWidth, innerHeight)
renderer.setPixelRatio(pr)
renderer.shadowMap.enabled = true
renderer.shadowMap.type = PCFSoftShadowMap


//////// LIGHT

const ambient = new HemisphereLight(
  new Color(1, 0.6, 0.2), new Color(0.2, 0.6, 1),
  1,
)
scene.add(ambient)

const light = new DirectionalLight(
  new Color(1, 0.9, 0.8),
  2,
)
light.position.set(-1, 1, 1)
scene.add(light)

light.castShadow = true
const { shadow } = light
shadow.mapSize.width  = 1024
shadow.mapSize.height = 1024
shadow.camera.near   =  0
shadow.camera.far    =  sqrt(3) * 2
shadow.camera.left   = -1
shadow.camera.right  =  1
shadow.camera.top    =  1
shadow.camera.bottom = -1
shadow.bias = -0.001
shadow.intensity = 0.62


//////// INTSTANCED CYLINDERS

const mat = new MeshLambertMaterial()
const geom = new CylinderGeometry(0.5, 0.5, 1, 24)
const obj = new InstancedMesh(geom, mat, N)
obj.scale.multiplyScalar(1.2)
obj.castShadow = true
obj.receiveShadow = true
obj.frustumCulled = false
obj.count = 0
scene.add(obj)

const rings = create_rings(64, 0.6, 1.2)
obj.add(rings)


const dummy = new Object3D()
const euler = new Euler()
const color = new Color()

const set_instance = (e: MessageEvent<number[]>) => {
  // max number of tests reached
  if (e.data.length === 2) {
    print(e.data[0], e.data[1])
    worker.terminate()
    return
  }

  const d = e.data,
     px = d[0], py = d[1], pz = d[2],
     rx = d[3], rz = d[4],
     sr = d[5], sh = d[6],
      n = d[7], tn = d[8]

  dummy.position.set(px, py, pz)
  dummy.setRotationFromEuler(euler.set(rx, 0, rz))
  dummy.scale.set(sr, sh, sr)
  dummy.updateMatrix()

  obj.setMatrixAt(obj.count, dummy.matrix)
  obj.setColorAt(obj.count, color
    .setRGB(
      rx / (TAU * 3) + sr / (MAX_R * 1.5),
      rz / (TAU * 2) + sh / (MAX_H * 2),
      1 - sr / (MAX_R * 3))
    .multiplyScalar(sh > (sr * 3) ? 2 : 1)
    .addScalar((sr / MAX_R + sh / MAX_H) < 0.38 ? 0.62 : 0)
  )

  obj.instanceMatrix.addUpdateRange(obj.count * 16, 16)
  obj.instanceColor!.addUpdateRange(obj.count *  3,  3)
  obj.instanceMatrix.needsUpdate = true
  obj.instanceColor!.needsUpdate = true
  obj.count = n

  print(n, tn)

  if (n === N) worker.terminate()
}


//////// PACKING WORKER

let worker: Worker

fetch("https://unpkg.com/simplex-noise@4.0.1/dist/esm/simplex-noise.js")
.then(res => res.text())
.then((noise: string) => {

  const worker_src = `
  const N = ${N}
  const MAX_R = ${MAX_R}
  const MIN_R = ${MIN_R}
  const MAX_H = ${MAX_H}
  const MIN_H = ${MIN_H}
  const MAX_TESTS = ${MAX_TESTS}
  ${noise.replace(/^export /gm, '')}
  ${worker_self.toString().replace(/^.+\n|\n.+$/g, "")}`

  worker = new Worker(URL.createObjectURL(new Blob(
    [ worker_src ],
    { type: "text/javascript" },
  )))

  // init packing
  worker.postMessage(optimized_cylinder())
  worker.onmessage = set_instance
})

function worker_self() {

  const {
    abs, cos, min, PI, random: rnd, sin, sqrt,
  } = Math
  const TAU = PI * 2
  const DIFF_R = MAX_R - MIN_R
  const DIFF_H = MAX_H - MIN_H

  // spatial hash grid dimensions
  const grid_side = 8
  const grid_area = grid_side * grid_side
  const cell_len = 1 / grid_side

  // init packing
  onmessage = (e) => {
    const geom = e.data
    const geom_len = geom.length

    const packed = new Float32Array(N * 7)
    const noise = createNoise2D()
    // packing matrix helpers
    const m  = new_mat4()
    const mp = new_mat4()


    // spatial hash grid
    const grid = new_grid()
    // hash processing helpers
    const cells_set = new Set
    const checked_ids = []

    // message values
    const msg = Array(9)
    let n = 0, tests_n = 0
    let
      // position x, y, z
      px = 0, py = 0, pz = 0,
      // scale radius, height
      sr = 0, sh = 0,
      // rotation x, z
      rx = 0, rz = 0

    // packing loop
    while (n < N) {
      sr = 0
      const n_ratio = 1 - n / N
      cells_set.clear()

      pack: while (sr === 0) {
        // position
        px = rnd_pos()
        py = rnd_pos()
        pz = rnd_pos()
        if (len2(px, pz) > 0.5 || abs(py) > 0.5) {
          sr = 0
          continue pack
        }

        // rotation
        const scl = 0.9
        rx = noise(py * scl,  pz * scl) * 0.5 + 0.5
        rz = noise(px * scl, -py * scl) * 0.5 + 0.5
        // smoothsteped
        const srx = smooth(rx)
        const srz = smooth(rz)
        rx *= TAU
        rz *= TAU

        // scale
        if (n < 100) {
          const n2 = n < 20 ? 2 : 1
          sr = MIN_R + srx * DIFF_R * n_ratio * n2
          sh = MIN_H + srz * DIFF_H * n_ratio * n2
        }
        else {
          sr = MIN_R + min(srx, rnd()) * DIFF_R * n_ratio
          sh = MIN_H + min(srz, rnd()) * DIFF_H * n_ratio
        }

        // compose packing matrix
        compose_mat4(m, px, py, pz, rx, rz, sr, sh)

        // apply matrix to a cloned basic geometry
        // check if the cylinder intersects a boundary
        const g = geom.slice()
        for (let i = 0; i < geom_len; i += 3) {
          apply_mat4(g, i, m)
          if (len2(g[i], g[i+2]) > 0.5 || abs(g[i+1]) > 0.5) {
            sr = 0
            cells_set.clear()
            continue pack
          }
          add_cell(g, i, cells_set)
        }

        checked_ids.length = 0
        // check meshes from hood cells
        for (const cell of cells_set) {
        const hood = grid[cell].hood
        for (const hood_cell of hood) {
        const packed_ids = grid[hood_cell].n
        for (const id of packed_ids) {

          if (tests_n === MAX_TESTS) {
            return postMessage([ n, tests_n ])
          }
          if (checked_ids.includes(id)) {
            continue
          }
          tests_n++
          checked_ids.push(id)

          const id7 = id * 7
          // compose packed matrix
          compose_mat4(
            mp,
            // px, py, pz
            packed[id7+0], packed[id7+1], packed[id7+2],
            // rx, rz
            packed[id7+3], packed[id7+4],
            // sr, sh
            packed[id7+5], packed[id7+6],
          )

          // apply packed matrix to a cloned geometry
          // of the basic cylinder
          const gp = geom.slice()
          for (let i = 0; i < geom_len; i += 3) {
            apply_mat4(gp, i, mp)
          }

          if (check_intersection(g, gp)) {
            sr = 0
            cells_set.clear()
            continue pack
          }
        }}}
      }

      // put cells ids of the cylinder into the grid hash
      for (const cell of cells_set!) {
        grid[cell].n.push(n)
      }

      const i = n * 7
      packed[i + 0] = msg[0] = px
      packed[i + 1] = msg[1] = py
      packed[i + 2] = msg[2] = pz
      packed[i + 3] = msg[3] = rx
      packed[i + 4] = msg[4] = rz
      packed[i + 5] = msg[5] = sr
      packed[i + 6] = msg[6] = sh
      msg[7] = ++n
      msg[8] = tests_n

      postMessage(msg)
    }
  }

  function check_intersection(g, gp) {
    // check if the curr cylinder intersects the packed
    // both are projected to 2d planes
    for (let i = 0; i < gp.length; i += 9) {
      const x1o = gp[i+0], x2o = gp[i+3], x3o = gp[i+6]
      const y1o = gp[i+1], y2o = gp[i+4], y3o = gp[i+7]
      const z1o = gp[i+2], z2o = gp[i+5], z3o = gp[i+8]

    for (let j = 0; j < g.length; j += 9) {
      const x1 = g[j+0], x2 = g[j+3], x3 = g[j+6]
      const y1 = g[j+1], y2 = g[j+4], y3 = g[j+7]
      const z1 = g[j+2], z2 = g[j+5], z3 = g[j+8]

    if (
      check_triangles(
        x1o, y1o, x2o, y2o, x3o, y3o,
        x1,  y1,  x2,  y2,  x3,  y3) ||
      check_triangles(
        y1o, z1o, y2o, z2o, y3o, z3o,
        y1,  z1,  y2,  z2,  y3,  z3) ||
      check_triangles(
        z1o, x1o, z2o, x2o, z3o, x3o,
        z1,  x1,  z2,  x2,  z3,  x3)
    ) {
      return true
    }}}
    return false
  }

  function smooth(t) {
    return t * t * (3 - 2 * t)
  }

  function rnd_pos() {
    return rnd() * 2 - 1
  }

  function len2(a, b) {
    return sqrt(a * a + b * b)
  }

  function new_mat4() {
    return new Float32Array([
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
  ])}

  function compose_mat4(
    m,
    pos_x, pos_y, pos_z,
    rot_x, rot_z,
    scl_r, scl_h,
  ) {
    const
      cx = cos(rot_x), sx = sin(rot_x),
      cz = cos(rot_z), sz = sin(rot_z),
      cxz = cx * cz, cxsz = cx * sz,
      sxz = sx * sz, sxcz = sx * cz
    m[0] =   cz * scl_r
    m[1] = cxsz * scl_r
    m[2] =  sxz * scl_r
    m[4] =  -sz * scl_h
    m[5] =  cxz * scl_h
    m[6] = sxcz * scl_h
    m[9] =  -sx * scl_r
    m[10] =  cx * scl_r
    m[12] = pos_x
    m[13] = pos_y
    m[14] = pos_z
  }

  function apply_mat4(g, i, m) {
    const x = g[i], y = g[i+1], z = g[i+2]
    g[i+0] = m[12] + m[0] * x + m[4] * y
    g[i+1] = m[13] + m[1] * x + m[5] * y + m[9]  * z
    g[i+2] = m[14] + m[2] * x + m[6] * y + m[10] * z
  }
  
  function check_triangles(
    x0, y0, x1, y1, x2, y2,
    x3, y3, x4, y4, x5, y5,
  ) {
    const arg = arguments
    // AABB check, early return
    let min1x = x0, max1x = min1x,
        min1y = y0, max1y = min1y,
        min2x = x3, max2x = min2x,
        min2y = y3, max2y = min2y
    for (let i = 2, j = 8; i < 6; i += 2, j += 2) {
      const p1x = arg[i], p1y = arg[i+1]
      const p2x = arg[j], p2y = arg[j+1]
      if (p1x < min1x) min1x = p1x
      if (p1x > max1x) max1x = p1x
      if (p1y < min1y) min1y = p1y
      if (p1y > max1y) max1y = p1y
      if (p2x < min2x) min2x = p2x
      if (p2x > max2x) max2x = p2x
      if (p2y < min2y) min2y = p2y
      if (p2y > max2y) max2y = p2y
    }
    if (
      max1x < min2x || max2x < min1x ||
      max1y < min2y || max2y < min1y
    ) {
      return false
    }
    // SAT on 6 axes
    const axes = [
      y0 - y1, x1 - x0, y1 - y2,
      x2 - x1, y2 - y0, x0 - x2,
      y3 - y4, x4 - x3, y4 - y5,
      x5 - x4, y5 - y3, x3 - x5,
    ]
    for (let i = 0; i < 12; i += 2) {
      const ax = axes[i], ay = axes[i+1]
      let min1d = x0 * ax + y0 * ay,
          min2d = x3 * ax + y3 * ay,
          max1d = min1d,
          max2d = min2d
      for (let i = 2, j = 8; i < 6; i += 2, j += 2) {
        const d1 = arg[i] * ax + arg[i+1] * ay
        const d2 = arg[j] * ax + arg[j+1] * ay
        if (d1 < min1d) min1d = d1;
        if (d2 < min2d) min2d = d2
        if (d1 > max1d) max1d = d1;
        if (d2 > max2d) max2d = d2
      }
      if (max1d < min2d || max2d < min1d) {
        return false
      }
    }
    return true
  }

  function add_cell(g, i, set) {
    const x = (g[i+0] + 0.5) / cell_len | 0,
          y = (g[i+1] + 0.5) / cell_len | 0,
          z = (g[i+2] + 0.5) / cell_len | 0,
          cell_id = y * grid_area + z * grid_side + x
    set.add(cell_id)
  }

  function new_grid() {
    const N = grid_area * grid_side
    const grid = Array(N)

    for (let i = 0; i < N; i++) {
      const x = i % grid_side
      const y = i / grid_area | 0
      const z = (i % grid_area) / grid_side | 0

      // array of ids of packed cylinders
      // related to current (x,y,z) cell
      const n = []
      // set of cells related to current (x,y,z) cell
      const hood = new Set
      set_hood(hood, x, y, z, { xb: grid_side })
      set_hood(hood, x, y, z, { yb: grid_side })
      set_hood(hood, x, y, z, { zb: grid_side })

      grid[i] = { hood, n }
    }

    return grid
  }

  function set_hood(
    hood, x, y, z,
    { xb = 1, yb = 1, zb = 1 },
  ) {
    for (let ix = -xb; ix <= xb; ix++)
    for (let iy = -yb; iy <= yb; iy++)
    for (let iz = -zb; iz <= zb; iz++) {
      const nx = x + ix, ny = y + iy, nz = z + iz
      if (
        nx < 0 || ny < 0 || nz < 0 ||
        nx >= grid_side ||
        ny >= grid_side ||
        nz >= grid_side
      ) continue
      hood.add(ny * grid_area + nz * grid_side + nx)
    }
  }
}


//////// ANIMATION LOOP

const DUR   = 3000
const DELAY = 2000

const get_next_qt = (() => {
  const new_qt = (v: Vector3, a: number) => (
    new Quaternion().setFromAxisAngle(v, a)
  )

  const HPI = PI / 2

  const vx = new Vector3(1, 0, 0)
  const qtxs = [
    new_qt(vx, -HPI),
    new_qt(vx, 0),
    new_qt(vx, HPI),
  ]

  const vy = new Vector3(0, 1, 0)
  const qtys = [
    new_qt(vy, -HPI),
    new_qt(vy, 0),
    new_qt(vy, HPI),
  ]

  // 0b0101 gives zero rotation, qtxs[1] * qtys[1]
  let prev_qt_id = 5
  let xqt: Quaternion
  let yqt: Quaternion

  const set_rnd_qt = () => {
    const id = ((rnd() * 3) << 2) | ((rnd() * 3) & 3)
    // don't return zero or same rotation
    if (id === prev_qt_id || id === 5) {
      return set_rnd_qt()
    }
    prev_qt_id = id
    xqt = qtxs[id >> 2]
    yqt = qtys[id & 3]
  }

  return (prev_qt: Quaternion) => {
    set_rnd_qt()
    return prev_qt.clone()
      .premultiply(xqt)
      .premultiply(yqt)
  }
})()

const qt = obj.quaternion

// quadinout
const ease = (t: number) => t < 0.5
  ? 2 * t * t
  : t * 4 - 2 * t * t - 1

const rnd_rot = () => {
  const prev_qt = qt.clone()
  const next_qt = get_next_qt(prev_qt)

  const start = performance.now()

  const tick = () => {
    const now = performance.now()
    const t = ease(min(1, (now - start) / DUR))
    qt.slerpQuaternions(prev_qt, next_qt, t)
    if (t < 1) requestAnimationFrame(tick)
    else setTimeout(rnd_rot, DELAY)
  }

  requestAnimationFrame(tick)
}

setTimeout(rnd_rot, DELAY)


//////// RENDER

renderer.setAnimationLoop(() => {
  renderer.render(scene, camera)
})


//////// UTILS

function str(num: number) {
  return num.toLocaleString()
}

function print(n: number, tests_n: number) {
  log.textContent = `${str(n)} cylinders\n${str(tests_n)} tests`
}

function set_camera_aspect(camera: OrthographicCamera) {
  if (innerWidth > innerHeight) {
    const asp = innerWidth / innerHeight
    camera.left   = -asp
    camera.right  =  asp
    camera.top    =  1
    camera.bottom = -1
  } else {
    const asp = innerHeight / innerWidth
    camera.left   = -1
    camera.right  =  1
    camera.top    =  asp
    camera.bottom = -asp
  }
  camera.updateProjectionMatrix()
}

onresize = () => {
  set_camera_aspect(camera)
  renderer.setSize(innerWidth, innerHeight)
}

function optimized_cylinder() {
  const s = 0.5
  const a_len = PI / 3
  const buf: number[] = []

  let x = cos(a_len) * s
  let z = sin(a_len) * s
  for (let i = -1; i < 2; i += 2) {
    const y = s * i
    buf.push(
      s,  y,  0,  x,  y,  z,  x,  y, -z,
     -s,  y,  0, -x,  y, -z, -x,  y,  z,
      x,  y,  z, -x,  y,  z,  x,  y, -z,
     -x,  y, -z,  x,  y, -z, -x,  y,  z,
   )
  }

  for (let i = 0; i < 3; i++) {
    const a = a_len * i
    x = cos(a) * s
    z = sin(a) * s
    buf.push(
      x,  s,  z, -x,  s, -z, -x, -s, -z,
     -x, -s, -z,  x, -s,  z,  x,  s,  z
    )
  }

  return buf
}

function create_rings(
  segments: number,
  radius: number,
  height: number,
) {
  const pos_buf = new Float32Array(segments * 12)
  const a_len = PI / segments * 2
  height *= 0.5

  for (let i = 0; i < segments; i++) {
    let j = i * 12
    for (let y = -1; y < 2; y += 2) {
      j += (y + 1) * 3
      for (let a = 0; a < 2; a++) {
        const ai = (i + a) * a_len
        const k = j + a * 3
        pos_buf[k+0] = cos(ai) * radius
        pos_buf[k+2] = sin(ai) * radius
        pos_buf[k+1] = y * height
      }
    }
  }

  const geom = new BufferGeometry()
  const attr = new BufferAttribute(pos_buf, 3)
  geom.setAttribute("position", attr)
  const mat = new LineBasicMaterial()
  return new LineSegments(geom, mat)
}

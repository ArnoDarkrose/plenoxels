use glam::{Vec3, Mat4, Vec4, UVec3};
use image::{ImageBuffer, Rgba};
use crate::Cell;

const SH_WIDTH: usize = 9;

pub fn launch(width: u32, height: u32, proj_inv_matrix: Mat4, world_view_inv: Mat4, grid: &[Cell]) -> ImageBuffer<Rgba<u8>, Vec<u8>> { 

    ImageBuffer::from_par_fn(width, height, |x, y| {
        let (ray_dir, ray_pos) = get_ray_dir_pos(
            (x as f32 + 0.5) / width as f32, (y as f32 + 0.5) / height as f32, proj_inv_matrix, world_view_inv
        );

        //TODO changed box_min and box_max from -1.0 and 1.0
        let (near, far) = ray_box_intersect(ray_pos, ray_dir, Vec3::splat(-63.0), Vec3::splat(64.0));

        let color = ray_march_sh(grid, near, far, 128, ray_dir);

        Rgba([(color.x * 255.0) as u8, (color.y * 255.0) as u8, (color.z * 255.0) as u8, 255])
    })
}

fn ray_march_sh(grid: &[Cell], tmin: f32, tmax: f32, grid_size: usize, ray_dir: Vec3) -> Vec3 {
    let dt = 1.0 / grid_size as f32;

    let mut t = tmin;

    //let mut transmitted_light = 1.0;
    //let mut steps = 0;

    //let mut transmitted_light = vec![1.0];
    let mut last_transmitted_light = 1.0;

    let mut color = Vec3::splat(0.0);
    
    while t < tmax {
        let cur_cell = grid_sample(grid, ray_dir * t, UVec3::splat(grid_size as u32));

        /*
        if cur_cell.density > 0.0{
            return Vec3::splat((0..=steps).map(|i| transmitted_light[i] * (1.0 - f32::exp(-cur_cell.density * dt))).sum::<f32>());
            return Vec3::splat(cur_cell.density);
        }
        */
        
        color += last_transmitted_light * (1.0 - f32::exp(-cur_cell.density * dt)) * Vec3::new(
            eval_sh(&cur_cell.r_sh, ray_dir),
            eval_sh(&cur_cell.g_sh, ray_dir),
            eval_sh(&cur_cell.b_sh, ray_dir)
        );

        last_transmitted_light *= f32::exp(-cur_cell.density * dt);
        //transmitted_light.push(last_transmitted_light);

        t += dt;
        //steps += 1;
    }

    color
}

fn sh_eval_2(ray_dir: Vec3) -> [f32; SH_WIDTH] {
    let x = ray_dir.x;
    let y = ray_dir.y;
    let z = ray_dir.z;

    let c0 = x;
    let s0 = y;

    let z2 = z*z;

    let mut res = [0.0; SH_WIDTH];
    res[0] = 0.28209479177387814;
    res[2] = z * 0.488602511902919923;
    res[6] = z2 * 0.94617469575756008 + -0.315391565252520045;

    let tmp_a = -0.488602511902919978;
    res[3] = tmp_a * c0;
    res[1] = tmp_a * s0;

    let tmp_b = z * -1.09254843059207896;
    res[7] = tmp_b * c0;
    res[5] = tmp_b * s0;

    let c1 = x * c0 - y * s0;
    let s1 = x * s0 + y * c0;

    let tmp_c = 0.546274215296039478;
    res[8] = tmp_c * c1;
    res[4] = tmp_c * s1;

    res
}
fn eval_sh(sh: &[f32; SH_WIDTH], ray_dir: Vec3) -> f32 { 
    sh_eval_2(ray_dir).into_iter().enumerate().map(|(i, v)| v * sh[i]).sum()
}

fn grid_sample(grid: &[Cell], position: Vec3, grid_dimensions: UVec3) -> Cell {
    let x = position.x.clamp(0.0, grid_dimensions.x as f32 - 1.0);
    let y = position.y.clamp(0.0, grid_dimensions.y as f32 - 1.0);
    let z = position.z.clamp(0.0, grid_dimensions.z as f32 - 1.0);

    let x0 = x.floor() as usize;
    let y0 = y.floor() as usize;
    let z0 = z.floor() as usize;
    let x1 = (x0 + 1).min(grid_dimensions.x as usize - 1);
    let y1 = (y0 + 1).min(grid_dimensions.y as usize - 1);
    let z1 = (z0 + 1).min(grid_dimensions.z as usize - 1);

    let xd = x - x0 as f32;
    let yd = y - y0 as f32;
    let zd = z - z0 as f32;

    let index = |x, y, z| x + grid_dimensions.x as usize * (y + grid_dimensions.y as usize * z);

    let c000 = &grid[index(x0, y0, z0)];
    let c100 = &grid[index(x1, y0, z0)];
    let c010 = &grid[index(x0, y1, z0)];
    let c110 = &grid[index(x1, y1, z0)];
    let c001 = &grid[index(x0, y0, z1)];
    let c101 = &grid[index(x1, y0, z1)];
    let c011 = &grid[index(x0, y1, z1)];
    let c111 = &grid[index(x1, y1, z1)];

    let interpolate = |c0: &Cell, c1: &Cell, t: f32| -> Cell {
        let mut result = Cell::default();
        result.density = c0.density * (1.0 - t) + c1.density * t;
        for i in 0..SH_WIDTH {
            result.r_sh[i] = c0.r_sh[i] * (1.0 - t) + c1.r_sh[i] * t;
            result.g_sh[i] = c0.g_sh[i] * (1.0 - t) + c1.g_sh[i] * t;
            result.b_sh[i] = c0.b_sh[i] * (1.0 - t) + c1.b_sh[i] * t;
        }
        result
    };

    let c00 = interpolate(c000, c100, xd);
    let c01 = interpolate(c001, c101, xd);
    let c10 = interpolate(c010, c110, xd);
    let c11 = interpolate(c011, c111, xd);

    let c0 = interpolate(&c00, &c10, yd);
    let c1 = interpolate(&c01, &c11, yd);

    let ans = interpolate(&c0, &c1, zd);

    //TODO applied relu to density
    Cell {density: ans.density.max(0.0), ..ans}
}
#[inline]
fn eye_ray_dir(x: f32, y: f32, proj_inv_matrix: Mat4) -> Vec3{
    let mut pos = Vec4::new(2.0 * x - 1.0, 2.0 * y - 1.0, 0.0, 1.0);

    pos = proj_inv_matrix * pos;
    pos /= pos.w;

    Vec3::new(pos.x, pos.y, pos.z).normalize()
}

#[inline]
fn get_ray_dir_pos(x: f32, y: f32, proj_inv_matrix: Mat4, world_view_inv: Mat4) -> (Vec3, Vec3) {
    let ray_dir = eye_ray_dir(x, y, proj_inv_matrix);
    let ray_pos = Vec3::splat(0.0);

    let ray_pos = world_view_inv * ray_pos.extend(1.0);
    let ray_dir = world_view_inv * ray_dir.extend(0.0).normalize();

    (ray_dir.truncate(), ray_pos.truncate())
}

fn ray_box_intersect(ray_pos: Vec3, mut ray_dir: Vec3, box_min: Vec3, box_max: Vec3) -> (f32, f32) {
    ray_dir.x = 1.0 / ray_dir.x;
    ray_dir.y = 1.0 / ray_dir.y;
    ray_dir.z = 1.0 / ray_dir.z;

    let lo = ray_dir.x * (box_min.x - ray_pos.x);
    let hi = ray_dir.x * (box_max.x - ray_pos.x);
    
    let mut tmin = lo.min(hi);
    let mut tmax = lo.max(hi);

    let lo = ray_dir.y * (box_min.y - ray_pos.y);
    let hi = ray_dir.y * (box_max.y - ray_pos.y);

    tmin = lo.min(hi).max(tmin);
    tmax = lo.max(hi).min(tmax);

    let lo = ray_dir.z * (box_min.z - ray_pos.z);
    let hi = ray_dir.z * (box_max.z - ray_pos.z);

    tmin = lo.min(hi).max(tmin);
    tmax = lo.max(hi).min(tmax);

    (tmin, tmax)
}

fn _ray_march_const_fog(tmin: f32, tmax: f32) -> (Vec4, f32) {
    let dt = 0.05; //i suppose dt is 1/grid_size, so grid_size here is 20
    let mut t = tmin;
    let mut alpha = 1.0;

    let mut color = Vec4::splat(0.0);

    while t < tmax && alpha > 0.01 {
        let a = 0.025;
        color += a * alpha * Vec4::new(1.0, 1.0, 0.0, 0.0);
        alpha *= 1.0 - a;
        t += dt;
    }

    (color, alpha)
}

/*
fn ray_march_sh(grid: &[Cell], camera_pos: Vec3, ray_dir: Vec3, grid_dimensions: UVec3, max_steps: usize, step_size: f32) -> Vec3 {
    let mut position = camera_pos;
    for _ in 0..max_steps {
        if position.x < 0.0 || position.y < 0.0 || position.z < 0.0 || position.x >= grid_dimensions.x as f32 || position.y >= grid_dimensions.y as f32 || position.z >= grid_dimensions.z as f32 {
            break; // Exit if the ray is outside the grid bounds
        }
        let cell = grid_sample(grid, position, grid_dimensions); // Assume grid_sample is defined elsewhere
        if cell.density > 0.0 {
            // Assuming eval_sh is defined elsewhere and computes the color based on spherical harmonics coefficients
            let color = Vec3::new(
                eval_sh(&cell.r_sh, ray_dir),
                eval_sh(&cell.g_sh, ray_dir),
                eval_sh(&cell.b_sh, ray_dir),
            );
            return color; // Return the color if a non-transparent cell is hit
        }
        position += ray_dir * step_size; // Move the ray forward
    }
    Vec3::ZERO // Return black if no hit
}
*/

/*
fn generate_image(grid: &[Cell], grid_dimensions: UVec3, camera_pos: Vec3, view_dir: Vec3, resolution: UVec2, fov: f32) -> Vec<Vec3> {
    let aspect_ratio = resolution.x as f32 / resolution.y as f32;
    let proj_matrix = Mat4::perspective_rh_gl(fov.to_radians(), aspect_ratio, 0.1, 100.0);
    let view_matrix = Mat4::look_at_rh(camera_pos, camera_pos + view_dir, Vec3::Y);

    let mut image = Vec::with_capacity((resolution.x * resolution.y) as usize);

    for y in 0..resolution.y {
        for x in 0..resolution.x {
            let ndc_x = (x as f32 + 0.5) / resolution.x as f32 * 2.0 - 1.0;
            let ndc_y = (y as f32 + 0.5) / resolution.y as f32 * 2.0 - 1.0;
            let ray_clip = Vec4::new(ndc_x, ndc_y, -1.0, 1.0);
            let ray_eye = proj_matrix.inverse() * ray_clip;
            let ray_world = (view_matrix.inverse() * Vec4::new(ray_eye.x, ray_eye.y, -1.0, 0.0)).truncate().normalize();

            let color = ray_march_sh(grid, camera_pos, ray_world, grid_dimensions, 1000, 0.1);
            image.push(color);
        }
    }

    image
}
*/

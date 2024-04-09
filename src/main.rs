pub mod cpu;
use glam::{Vec3, Mat4, Vec4, UVec2};
use std::{
    sync::Arc,
    io::{BufReader, Read},
    fs::File,
    cell::RefCell,
    ptr
};

const SH_WIDTH: usize = 9;

#[derive(Debug, Clone, PartialEq)]
pub struct Cell {
    pub density: f32,
    pub r_sh: [f32; SH_WIDTH],
    pub g_sh: [f32; SH_WIDTH],
    pub b_sh: [f32; SH_WIDTH]
}

impl Default for Cell {
    fn default() -> Self {
        Cell {
            density: 0.0,
            r_sh: [0.0; 9],
            g_sh: [0.0; 9],
            b_sh: [0.0; 9]
        }
    }
}

impl Cell {
    pub fn new_safe(reader: Arc<RefCell<BufReader<File>>>) -> Self {
        let mut reader = reader.borrow_mut();

        let mut buf = [0; 4];

        reader.read_exact(&mut buf).expect("failed to read cell");
        let density = f32::from_le_bytes(buf);

        let mut r_sh = [0.0; SH_WIDTH];
        let mut g_sh = [0.0; SH_WIDTH];
        let mut b_sh = [0.0; SH_WIDTH];

        for i in 0..SH_WIDTH {
            reader.read_exact(&mut buf).expect("failed to read cell");
            r_sh[i] = f32::from_le_bytes(buf);
        }

        for i in 0..SH_WIDTH {
            reader.read_exact(&mut buf).expect("failed to read cell");
            g_sh[i] = f32::from_le_bytes(buf);
        }

        for i in 0..SH_WIDTH {
            reader.read_exact(&mut buf).expect("failed to read cell");
            b_sh[i] = f32::from_le_bytes(buf);
        }

        Self {density, r_sh, g_sh, b_sh}
    }

    pub fn new(reader: Arc<RefCell<File>>) -> Self {
        let mut reader = reader.borrow_mut();

        let mut buf = [0; 112];

        reader.read_exact(&mut  buf).expect("failed to read cell");


        let density;
        let r_sh;
        let g_sh;
        let b_sh;

        let buf:*const f32 = ptr::addr_of!(buf).cast();

        //density is placed at the same place as that the buf ponts at. r_sh, g_sh, and b_sh are places after that, so
        //i skeep 1, 1 + SH_WIDTH, 1 + 2 * SH_WIDTH, respectively to get them
        unsafe {
            density = *buf;

            r_sh = *buf.add(1).cast();

            g_sh = *buf.add(1 + SH_WIDTH).cast();

            b_sh = *buf.add(1 + SH_WIDTH*2).cast();
        }
       
        Self {density, r_sh, g_sh, b_sh}
    }
}

fn main() {
    let resolution = UVec2::new(512, 512);

    let proj_matrix = Mat4::perspective_rh(90.0_f32.to_radians(), 1.0, 0.1, 100.0).inverse();

    let mut grid = Vec::with_capacity(128_usize.pow(3));
    //let reader = Arc::new(RefCell::new(BufReader::new(File::open("model.dat").expect("failed to open file"))));

    //TODO i don't need Arc<RefCell> here
    let reader = Arc::new(RefCell::new(File::open("model.dat").expect("failed to open model file")));

    for _ in 0..128_usize.pow(3) {
        grid.push(Cell::new(Arc::clone(&reader)));
    }

    //TODO change 10 to 360
    for angle_y in (0..10).step_by(10).map(|v| v as f32) {
        let rotate_matrix = Mat4::from_rotation_y(angle_y.to_radians());

        let cam_pos = rotate_matrix * Vec4::new(0.0, 0.0, -3.0, 0.0) + Vec4::new(0.0, 1.5, 0.0, 1.0);

        let view = Mat4::look_at_rh(
            cam_pos.truncate(), Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)
        ).inverse();

        let img_buf = cpu::ray_march::launch(resolution.x, resolution.y, proj_matrix, view, &grid[..]);    
        img_buf.save(format!("images/{angle_y} degrees.png")).unwrap();
        
        println!("{angle_y}");
        /*
        let img_buf: Vec<_> = generate_image(
            &grid[..], 
            UVec3::new(128, 128, 128), 
            cam_pos.truncate(), 
            Vec3::new(0.0, 0.0, 0.0), 
            resolution, 90.0
        ).into_iter().map(|c| [(c.x * 255.0) as u8, (c.y * 255.0) as u8, (c.z * 255.0) as u8, 255]).flatten().collect();

        let image = RgbaImage::from_raw(resolution.x, resolution.y, img_buf).unwrap();

        image.save(format!("images/{angle_y} degrees.png")).unwrap();
        */
    }
}       

//showed that the file stores exactly 2097152 cells. This is 128^3
#[cfg(test)]
mod tests {
    use std::io::Read;
    #[test]
    fn test() -> std::io::Result<()> {
        let reader = std::io::BufReader::new(std::fs::File::open("model.dat")?);

        let mut buf = [0; 112];

        let shared = std::sync::Arc::new(std::cell::RefCell::new(reader));
        
        let mut counter = 0;

        while let Ok(_) = shared.borrow_mut().read_exact(&mut buf) {
            counter += 1;
        }

        println!("{counter}");

        Ok(())
    }

    use crate::Cell;
    use std::sync::Arc;
    #[test]
    #[ignore]
    fn test_read() -> std::io::Result<()>{
        let start = std::time::Instant::now();
        let mut grid = Vec::with_capacity(128_usize.pow(3));
        let reader = Arc::new(std::cell::RefCell::new(std::io::BufReader::new(std::fs::File::open("model.dat")?)));
        
        for _ in 0..128_usize.pow(3) {
            grid.push(Cell::new_safe(Arc::clone(&reader)));
        }

        let took_time = std::time::Instant::now().duration_since(start);
        dbg!(took_time);

        println!("{}", grid.len());

        println!("{:?}", grid[0]);
        Ok(())
    }

    #[test]
    #[ignore]
    fn test_unsafe_read() -> std::io::Result<()> {
        let start = std::time::Instant::now();
        let mut grid = Vec::with_capacity(128_usize.pow(3));
        let reader = Arc::new(std::cell::RefCell::new(std::fs::File::open("model.dat")?));
        
        for _ in 0..128_usize.pow(3) {
            grid.push(Cell::new(Arc::clone(&reader)));
        }

        let took_time = std::time::Instant::now().duration_since(start);
        dbg!(took_time);

        println!("{}", grid.len());
        println!("{:?}", grid[0]);
        Ok(())
    }

}